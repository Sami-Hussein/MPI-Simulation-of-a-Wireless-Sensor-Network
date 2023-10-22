
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <unistd.h>
#include <time.h>
#include <mpi.h>
#include "wsn.h"
// function to update the DateTime structure to represent 10 seconds time passage
void update_datetime(struct DateTime *dt) {
    dt->second += 10;
    dt->minute += dt->second / 60;
    dt->second %= 60;
    dt->hour += dt->minute / 60;
    dt->minute %= 60;
    
    // Increment the date components if necessary
    if (dt->hour >= 24) {
        dt->hour -= 24;
        dt->day++;
    }
}
void master(MPI_Comm master_comm, MPI_Comm slave_comm, int my_rank, int size, int nrows, int ncols, MPI_Datatype mpi_report, MPI_Datatype mpi_summary, int runtime);
void slave(MPI_Comm master_comm, MPI_Comm slave_comm, int my_rank, int base_station, MPI_Datatype mpi_report, MPI_Datatype mpi_summary);

int main(int argc, char *argv[])
{   
    // MPI setup
    int ndims = 2, size, my_rank, reorder, ierr;
    int nrows, ncols;
    int base_station;
    MPI_Comm comm2D;
    int dims[ndims];
    int wrap_around[ndims];
    int provided;
    // MPI_Init_threads for hybrid parallelism
    MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, &provided);
    
    // create datetime MPI_object
    MPI_Datatype mpi_datetime;
    MPI_Type_contiguous(6, MPI_INT, &mpi_datetime);
    MPI_Type_commit(&mpi_datetime);

    // Create MPI datatype for Report
    MPI_Datatype mpi_report;
    int block_lengths[8] = {1, 1, 1, 1, 1, 1, 1, 4}; // Added 4 for the distant_nbrs array
    MPI_Aint offsets[8];
    MPI_Datatype types[8] = {mpi_datetime, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_INT}; 
    offsets[0] = 0;
    offsets[1] = offsetof(struct Report, messages_counter);
    offsets[2] = offsetof(struct Report, nbrs_count);
    offsets[3] = offsetof(struct Report, reporting_node);
    offsets[4] = offsetof(struct Report, iteration);
    offsets[5] = offsetof(struct Report, time_taken);
    offsets[6] = offsetof(struct Report, start_time);
    offsets[7] = offsetof(struct Report, nbrs);
    MPI_Type_create_struct(8, block_lengths, offsets, types, &mpi_report);
    MPI_Type_commit(&mpi_report);

    // create MPI datatype for summary
    MPI_Datatype mpi_summary;
    int summary_lengths[2] = {1,1};
    MPI_Aint summary_offsets[2];
    MPI_Datatype summary_types[2] = {MPI_INT, MPI_INT};
    summary_offsets[0] = 0;
    summary_offsets[1] = offsetof(struct Summary, total_master_msgs);
    MPI_Type_create_struct(2, summary_lengths, summary_offsets, summary_types, &mpi_summary);
    MPI_Type_commit(&mpi_summary);

    // get comm size and current rank
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    base_station = size - 1; // base station is the last process
    int runtime = 0;            // simulation runtime decided by user
    
    /* process command line arguments*/
	if (argc == 4) {
		nrows = atoi (argv[1]);
		ncols = atoi (argv[2]);
        runtime = atoi(argv[3]);
		dims[0] = nrows; 
		dims[1] = ncols; 
		if( (nrows*ncols) + 1 != size) {
			if( my_rank ==0 ) printf("ERROR: nrows*ncols + 1)=%d * %d + 1 = %d != %d\n", nrows, ncols, nrows*ncols,size);
			MPI_Finalize(); 
			return 0;
		}
	} 
    else 
    {
        printf("ERROR: Program args must be: nrows, ncols, simulation time(s)");
		MPI_Finalize();
        return 0;
	}

    // create the wrap around for mapping the 2d cartesian grid [false, false]
    wrap_around[0] = wrap_around[1] = 0; 
    reorder = 1;
    ierr =0;
    MPI_Comm_split(MPI_COMM_WORLD, my_rank==base_station , 0, &comm2D);

    // if base station then run master func, otherwise run slave
    if(my_rank == base_station){
        master(MPI_COMM_WORLD, comm2D, my_rank, size, nrows, ncols, mpi_report, mpi_summary, runtime);
    }
    else{
        MPI_Comm slave_comm;
        ierr = MPI_Cart_create(comm2D, ndims, dims, wrap_around, reorder, &slave_comm);
        if(ierr != 0) printf("ERROR[%d] creating CART\n",ierr);
        slave(MPI_COMM_WORLD, slave_comm, my_rank, base_station, mpi_report, mpi_summary);
    }        
    MPI_Type_free(&mpi_datetime);
    MPI_Type_free(&mpi_report);
    MPI_Type_free(&mpi_summary);
    MPI_Finalize();
    return 0;
}

/*
    The base station/master, functions as a single MPI process,monitors incoming reports from charging 
    nodes. It uses two threads. One for receiving reports and another sending a termination signal 
    and summarizing the simulation.Upon report reception, it checks neighboring nodes and suggests 
    a node if no reports arrive within a designated timeframe, or notifies nodes if no nearby 
    charging options are available. It logs key metrics, including simulation time, report count, 
    message sends, and contents. 

    @param master_comm MPI communication context for the master node
    @param slave_comm MPI communication context for the slave nodes
    @param my_rank Rank or identifier of the master node in the MPI communicator
    @param size Total number of nodes in the MPI communicator
    @param nrows Number of rows in the grid
    @param ncols Number of columns in the grid
    @param mpi_report MPI data type for the report received from the slave nodes
    @param mpi_summary MPI data type for the summary from the slave nodes
    @param runtime Total runtime duration of the simulation
    */
void master(MPI_Comm master_comm, MPI_Comm slave_comm, int my_rank, int size, int nrows, int ncols, MPI_Datatype mpi_report, MPI_Datatype mpi_summary, int runtime)
{   
    // use MPI_Wtime at the start of the process, time calculations between slave and master will be done relative to this point
    double start_time = MPI_Wtime();
    int i;
    // create logfile
    char buf[256] = {0};
    FILE *pFile;
    printf("Base Station Rank: %d. Comm Size: %d: Grid Dimension =[%d x %d] \n",my_rank, size, nrows, ncols);
    sprintf(buf, "base_station.txt");
    pFile = fopen(buf, "w");
    // delcare time calculation variables
    double total_base_comm_time = 0;        // total time spent in communicating with base station
    double total_nbrs_comm_time = 0;        // total time spent in communicating with slave nodes
    MPI_Request exit_request[size-1];       // MPI exit request
    MPI_Status exit_status[size-1];         // MPI exit status

    // declare the exit flag
    int exit_flag = 0;
    // create 2 threads, one for receving and logging, other is for sending the termination signal and summary writing
    #pragma omp parallel sections
    {
        #pragma omp section
        {   
            double master_comm_time;
            // create array of report struct
            struct Report report_array[MAX_REPORTS];
            int offset = 0;         // offset inside the report array

            // run while loop until termination message
            while(!exit_flag)
            {   
                // in each iteration we must probe to check whether to receive a report or not
                MPI_Status probe_status;
                int report_flag = 0;
                //MPI_Iprobe( int source , int tag , MPI_Comm comm , int* flag , MPI_Status* status);
                MPI_Iprobe(MPI_ANY_SOURCE, 0, master_comm, &report_flag, &probe_status);
                if(report_flag)
                {
                    // if we must receive a report then receive it using MPI_Recv
                    MPI_Recv( &report_array[offset],1 , mpi_report , MPI_ANY_SOURCE , 0, master_comm, MPI_STATUS_IGNORE);
                    
                    // Record the receive time and measure the time taken to send to the master node
                    master_comm_time = fabs((MPI_Wtime() - start_time) - report_array[offset].start_time);
                    // print the iteration, log time, and node data
                    fprintf(pFile, "------------------------------------------------------------------------------------\n");
                    fprintf(pFile, "Iteration: %d \n",report_array[offset].iteration);
                    fprintf(pFile, "Logged Time : %d-%d-%d %d:%d:%d\n", 
                            report_array[offset].datetime.year,
                            report_array[offset].datetime.month, report_array[offset].datetime.day,
                            report_array[offset].datetime.hour, report_array[offset].datetime.minute,
                            report_array[offset].datetime.second);
                            
                    fprintf(pFile, "Number of adjacent nodes: %d\n", report_array[offset].nbrs_count);
                    fprintf(pFile, "Availablity to be considered full: 0\n");
                    fprintf(pFile, "\nReporting Node\t\tCoord\t\tPort Value\t\tAvailable Ports\n");
                    int node_rank = report_array[offset].reporting_node;
                    int reporting_coords[2] = {node_rank/ncols, node_rank % ncols};
                    int max_ports = NUM_PORTS;
                    fprintf(pFile,"\t%d\t\t\t\t(%d,%d)\t\t\t%d\t\t\t\t 0\n", report_array[offset].reporting_node,reporting_coords[0], reporting_coords[1], max_ports);
                    fprintf(pFile, "\nAdjacent Node\t\tCoord\t\tPort Value\t\tAvailable Ports\n");

                    // print all of the adjacent nodes (neighbours of the reporting node)
                    for(i=0;i<4;i++)
                    {
                        if(report_array[offset].nbrs[i]>=0)
                        {
                            int reporting_rank = report_array[offset].nbrs[i];
                            int r = reporting_rank/ncols;
                            int c = reporting_rank%ncols;
                            fprintf(pFile, "\t%d\t\t\t\t(%d,%d)\t\t\t%d\t\t\t\t 0\n", reporting_rank, r, c,max_ports );
                        }
                    }
                    // print all of the nearby nodes
                    int nodes[8];
                    fprintf(pFile, "\nNearby Nodes\t\tCoord\n");

                    // calculate the ranks of all the nearby nodes and print to logfile
                    int dr[] = {-1, -1,  1, 1,  0, 0, -2, 2};
                    int dc[] = {-1,  1, -1, 1, -2, 2,  0, 0};
                    for (int i = 0; i < 8; i++) {
                        int dx = 0, dy = 0;
                        int rank;
                        dx = reporting_coords[0] -dr[i];
                        dy = reporting_coords[1] -dc[i];
                        rank = (dx * ncols) + dy;
                        nodes[i] = rank;

                        if(rank >= my_rank || rank < 0 || dx < 0 || dy < 0 || dx >= nrows || dy >= ncols)
                        {
                            rank = -1;
                        }
                        nodes[i] = rank;
                        if(rank != -1)
                            fprintf(pFile, "\t%d\t\t\t\t(%d,%d)\n",rank,dx,dy);
                    }

                    // check which nodes are available ie. did not submit a report in the last 2 iterations
                    int node_found = 0;         // flag to check whether we found at least 1 available node
                    int iteration_diff = 2;     // number of iteration difference between current index in this case its 2
                    
                    // iterate through all ranks calculated earlier
                    for(int i = 0; i< 8 && !node_found; i++)
                    {
                        // for each rank check whether it submitted a report in the last (iteration_diff) iterations
                        int distant_rank = nodes[i];
                        if(distant_rank >= 0)
                        {
                            for (int j = 0; j < offset; j++) {
                                int iter_diff = abs(report_array[j].iteration - report_array[offset].iteration);
                                if (report_array[j].reporting_node == distant_rank && iter_diff < iteration_diff)
                                {
                                    nodes[i] = -1;
                                }
                            }
                            
                        }
                    }
                    // print the available nodes and communication data
                    fprintf(pFile, "\nAvailable station nearby (no report received in the last %d iterations): ", iteration_diff);

                    for(i = 0; i < 8; i++){ 
                        if (nodes[i] != -1)
                            fprintf(pFile, "%d ", nodes[i]);
                    }
                    fprintf(pFile, "\nTotal Messages sent: %d\n", report_array[offset].messages_counter);
                    fprintf(pFile,"Communication Time With neighbours (s) : %f\n", report_array[offset].time_taken);
                    fprintf(pFile,"Communication Time With Master (s): %f\n",master_comm_time);
                    total_nbrs_comm_time += report_array[offset].time_taken;
                    total_base_comm_time += master_comm_time;

                    // Send the list of available nearby nodes to the reporting node as a suggestion
                    MPI_Request send_suggestion;
                    MPI_Isend(nodes , 8, MPI_INT , report_array[offset].reporting_node , 0 , master_comm , &send_suggestion);

                    // if last iteration then wait for send otherwise a process might be left hanging
                    if(exit_flag)
                        MPI_Wait( &send_suggestion , MPI_STATUS_IGNORE);
                    offset++;
                    fprintf(pFile, "\n"); // Add a separator between reports

                }
            }
        }

            #pragma omp section
            {   
                // sleep for t amount of seconds (simulation time) befor sending exit request
                sleep(runtime);
                int exit_all = EXIT_TAG;
                exit_flag = 1;
                // send exit message to all nodes
                for (int dest = 0; dest < my_rank; dest++)
                {
                    MPI_Isend(&exit_all, 1, MPI_INT, dest, EXIT_TAG, master_comm, &exit_request[dest]);
                }
                // wait for all messages to be sent
                MPI_Waitall( my_rank, exit_request, exit_status);

                printf("Ending Simulation, Creating Summary...\n");
                // create summary
                int total_master_calls = 0;         // number of calls to master node (total)
                int total_nbr_calls = 0;            // number of calls to neighbours (total)

                // receive a summary report from each slave node (for correctness checking and analysis)
                for(int i = 0; i < my_rank; i++)
                {
                    struct Summary summary;
                    MPI_Recv( &summary, 1 , mpi_summary , i , SUMMARY_TAG , master_comm , MPI_STATUS_IGNORE);
                    total_master_calls += summary.total_master_msgs;
                    total_nbr_calls += summary.total_nbr_msgs;
                }

                // print the summary to the logfile
                fprintf(pFile, "------------------------------------------------------------------------------------\n");
                fprintf(pFile, "\t\t\t\t\t\t\t\t\tSUMMARY\n");
                fprintf(pFile, "------------------------------------------------------------------------------------\n");

                fprintf(pFile, "Total Calls to Base Station: %d \n",total_master_calls);
                fprintf(pFile, "Total Calls to neighbours: %d \n",total_nbr_calls);
                fprintf(pFile, "Total Communication time with Base Station: %f \n",total_base_comm_time);
                fprintf(pFile, "Total Communication time with neighbours: %f \n",total_nbrs_comm_time);
                fprintf(pFile, "Average Communication time with Base Station: %f \n",total_base_comm_time/total_master_calls);
                fprintf(pFile, "Average Communication time with neighbours: %f \n",total_nbrs_comm_time/total_nbr_calls);
                fprintf(pFile, "Total Simulation time: %ds\n", runtime);
                fclose(pFile);
            }

    }
    

}

/*
    Charging station/Slave nodes function is used to represent the functionality of a single slave node in the WSN.
    Each node creates NUM_PORTS amount of threads, each with their own while loop to simulate time passage.
    Each iteration is followed by a 0.5s sleep/wait which represents 10s in simulated time.
    When a node has all its ports filled then it will query its neighbours for data and possibly notify the
    base station if all neighbours are full. Finally, the node submits a summary report for further analysis
    of the simulation.

    @param master_comm MPI communication context with the master node
    @param slave_comm MPI communication context specific to the slave nodes
    @param my_rank Rank or identifier of the current slave node in the MPI communicator
    @param base_station Identifier or address of the base station or master node
    @param mpi_report MPI data type for the report to be communicated to the master node
    @param mpi_summary MPI data type for the summary exchanged with the master node
    */
void slave(MPI_Comm master_comm, MPI_Comm slave_comm, int my_rank, int base_station, MPI_Datatype mpi_report, MPI_Datatype mpi_summary)
{
    // use MPI_Wtime at the start of each process to make all the time calculations with respect to it
    double node_start_time = MPI_Wtime();
    // time measurement variables declared
    double nbrs_start_time;             // measuring the start time when asking neighbours for data
    double nbrs_time_taken;             // measuring the total time when asking and receiving neighbour data
    // measure current time, the simulated time will be based on the current time
    time_t currentTime;
    struct tm *localTime;
    time(&currentTime);                     // Get current time in seconds since epoch
    localTime = localtime(&currentTime);    // Convert to local time

    // populate the datetime struct with the measured current time
    struct DateTime current_datetime = {
        .year = localTime->tm_year + 1900, 
        .month = localTime->tm_mon + 1,   
        .day = localTime->tm_mday,
        .hour = localTime->tm_hour,
        .minute = localTime->tm_min,
        .second = localTime->tm_sec
    };

    // create a cycle log array, whith each entry representing a log of the cycle information
    struct Cycle log_array[LOG_SIZE];
    int neighbours[4];                      // Neighbouring ranks: [left, right, top, bottom]
    int recv_data[4] = {-1,-1,-1,-1};       // data received from neighbours
    int ndims = 2;                          // number of dimensions
    int my_cart_rank;                       // node's rank in the cart
    int i;  
    // file buffer for debugging
    char buf[256] = {0};
    FILE *pFile;

    // Declaring MPI Requests and Statuses
    MPI_Request ask_request[4];             // Ask request used when asking neighbours for data
    MPI_Request respond_request;            // Request used by each neighbour when responding to the requesting node
    MPI_Status respond_status;              // Status of the response sent
    MPI_Request receive_request;            // receive request, used for receiving the neighbouring data
    MPI_Request share_data[4];              // share_data request, used to share a node's data before asking
    MPI_Request stop_req;                   // stop_req, used for receiving the termination msg from master node

    // Getting Cartesian grid details
    int coord[2];
    MPI_Cart_coords(slave_comm, my_rank, ndims, coord);                                 // get local coords
    MPI_Cart_rank(slave_comm, coord, &my_cart_rank);                                    // get current rank
    MPI_Cart_shift( slave_comm, SHIFT_ROW, DISP, &neighbours[2], &neighbours[3]);       // shift row
    MPI_Cart_shift( slave_comm, SHIFT_COL, DISP, &neighbours[0], &neighbours[1]);       // shift col

    // create a logfile for debugging 
    sprintf(buf, "log_%d.txt", my_rank);
    pFile = fopen(buf, "w");
    fprintf(pFile, "Global rank: %d. Cart rank: %d. Coord: (%d, %d). Left: %d. Right: %d. Top: %d. Bottom: %d\n\n", 
    my_rank, my_cart_rank, coord[0],coord[1], neighbours[0], neighbours[1], neighbours[2], neighbours[3]);

    // exit all flag for master termination msg
    int exit_all = 0;
    int master_message_count = 0;   // number of master calls (total)
    int neighbour_message_count = 0;    // number of neighbour calls (total)
    // receive termination message
    MPI_Irecv(&exit_all, 1, MPI_INT, base_station, EXIT_TAG, master_comm, &stop_req);
    // delcare pod_status and index var
    int index = 0;
    int pod_status = 0;
    // create NUM_PORTS amount of threads to simulate each thread
    #pragma omp parallel private(pod_status) num_threads(NUM_PORTS)
    {   
        // get thread id and seed the random numbers
        int thread_id = omp_get_thread_num();
        srand(time(NULL) + my_rank + thread_id);

        // each thread runs a while loop that gets syncrhonized during each iteration
        while(!exit_all) 
        {   
            // pod status represents whether a port/thread is in use or not 0 == not in use and 1 == in use
            pod_status = rand() % 2;            

            // when initially filling up the queue we let one thread create a struct and insert it
            if(index < LOG_SIZE){
                #pragma omp single
                {
                    struct Cycle mylog;
                    mylog.counter = 0;
                    mylog.cycle_num = index;
                    mylog.datetime = current_datetime;
                    log_array[index % LOG_SIZE] = mylog;
                }
            }

            // all threads must wait for the one thread to create the struct otherwise we cant write to it
            #pragma omp barrier

            // only one thread needs to update the time (each iteration is 10 seconds) and reset the port counter
            #pragma omp single
            {
                update_datetime(&current_datetime);
                if(index >= LOG_SIZE)
                {
                    // if we maxed out the circular queue then we need to overwrite the existing struct content
                    if(log_array[index % LOG_SIZE].cycle_num != index)
                    {   
                        // reset the counter, we effectively overwrite the oldest record
                        log_array[index % LOG_SIZE].counter = 0;
                    }          
                }
                // time is updated regardless
                log_array[index % LOG_SIZE].datetime = current_datetime;
            }

            // wait for one  thread to update time and reset counter (if needed)
            #pragma omp barrier                   
            log_array[index % LOG_SIZE].cycle_num = index;  // update the log index

            #pragma omp atomic
            log_array[index % LOG_SIZE].counter += pod_status;  // update the counter representing the number of ports in use

            #pragma omp barrier

                // one thread logs debugging information
                #pragma omp single
                {
                    // print log into the logfile
                    fprintf(pFile, "Log Entry %d:\nDate: %04d-%02d-%02d\nTime: %02d:%02d:%02d\nCounter: %u\n\n",
                    log_array[index % LOG_SIZE].cycle_num,
                    log_array[index % LOG_SIZE].datetime.year,
                    log_array[index % LOG_SIZE].datetime.month,
                    log_array[index % LOG_SIZE].datetime.day,
                    log_array[index % LOG_SIZE].datetime.hour,
                    log_array[index % LOG_SIZE].datetime.minute,
                    log_array[index % LOG_SIZE].datetime.second,
                    log_array[index % LOG_SIZE].counter);
                }
                // if all ports are in use then we need to ask for data
                if(log_array[index % LOG_SIZE].counter >= NUM_PORTS)
                {
                    
                    #pragma omp single 
                    {
                        int alerted_nbrs = 0;       // number of alerted neighbours
                        int total_nbrs = 0;         // total neighbours nearby
                        printf("P%d has no vacant ports.\n", my_rank);
                        fflush(stdout);

                        nbrs_start_time = MPI_Wtime();  // start measuring time for neighbour comms
                        
                        // ask all of the node's neighbours for data and send the current node's data to avoid a deadlock
                        for(i = 0; i < 4; i++)
                        {   
                            if(neighbours[i] >= 0){
                                total_nbrs++;
                                MPI_Isend(&index, 1, MPI_INT, neighbours[i], ASK_TAG, slave_comm, &ask_request[i]);
                                MPI_Isend(&log_array[index % LOG_SIZE].counter, 1, MPI_INT, neighbours[i], DATA_TAG, slave_comm, &share_data[i]);
                            }
                        }

                        // Check if neighbors exists and receive data if available
                        for (i = 0; i < 4; i++) {
                            if (neighbours[i] >= 0) {
                                    // A message is available, receive it
                                    MPI_Recv(&recv_data[i], 1, MPI_INT, neighbours[i], DATA_TAG, slave_comm, MPI_STATUS_IGNORE);
                                    alerted_nbrs++;         // increment the counter for the alerted neighbours
                            }
                        }

                        // iterate through the received data and check if all neighbours are occupied or not
                        int found_vacant = 0;
                        for(i = 0; i< 4; i++)
                        {
                            if(recv_data[i] >= 0 && recv_data[i] < 4){
                                recv_data[i] = 4 - recv_data[i];
                                printf("P%d suggests to use P%d with %d vacant ports\n", my_rank, neighbours[i], recv_data[i]);
                                fflush(stdout);
                                // if we found a match for a valid node then we just end at the first correct result
                                found_vacant = 1;
                            }
                        }

                        nbrs_time_taken = MPI_Wtime() - nbrs_start_time;

                        if(!found_vacant)
                        {   
                            printf("P%d Submits Report to Base Station.\n",my_rank);
                            // increment master_message counter
                            master_message_count++;
                            // increment neighbour counter for summary
                            neighbour_message_count += total_nbrs;
                            // create report
                            struct Report report;
                            report.datetime = log_array[index % LOG_SIZE].datetime;
                            report.iteration = index;                                       // iteration at which the error occured
                            report.reporting_node = my_rank;                                // rank of the reporting node
                            report.nbrs_count = total_nbrs;                                 // total number of neighbours for a node
                            report.time_taken = nbrs_time_taken;                            // time taken in neighbour communication
                            report.messages_counter = master_message_count;                 // number of calls made to master
                            report.start_time = (MPI_Wtime() - node_start_time); // time at which the communication with master started

                            for (int i = 0; i < 4; i++) {
                                report.nbrs[i] = neighbours[i];
                            }
                            // send report to base station
                            MPI_Send( &report , 1 , mpi_report , base_station , 0, master_comm);

                            // create int array to receive the suggested nodes
                            int suggested_nodes[8];
                            MPI_Recv( suggested_nodes , 8 , MPI_INT, base_station , 0, master_comm , MPI_STATUS_IGNORE);

                            // print the suggested nodes to terminal
                            printf("P%d Base Station suggests: ", my_rank);
                            fflush(stdout);
                            int sugg_count = 0;
                            for(i=0;i<8;i++)
                            {
                                if(suggested_nodes[i] != -1)    // check whether it is a valid node
                                {   sugg_count++;
                                    printf("%d ", suggested_nodes[i]);
                                    fflush(stdout);
                                }
                            }
                            if(sugg_count == 0) // if no node has been suggested then indicate in terminal
                                {
                                    printf("No Nearby Nodes Available.");
                                    fflush(stdout);
                                }
                            printf("\n" );
                            fflush(stdout);
                        }
                    }
                }
            
            #pragma omp barrier
            usleep(500000); // sleep for 0.5s to simulate 10s in the simulated time
            #pragma omp single
            {
                // test end condition
                MPI_Test( &stop_req, &exit_all , MPI_STATUS_IGNORE);
                MPI_Status probe_status;
                // testing for end condition must be done before sending to choose blocking/nonblocking send 

                // probe to check whether a ask/request has been received from a neighbour
                int send_flag = 0;
                int recv_index;
                MPI_Iprobe(MPI_ANY_SOURCE, ASK_TAG, slave_comm, &send_flag, &probe_status);
                if(send_flag)
                {
                    // if a neighbouring node asked for data, then we must send our data
                    MPI_Irecv(&recv_index, 1 , MPI_INT , probe_status.MPI_SOURCE , ASK_TAG , slave_comm, &receive_request );
                    MPI_Isend(&log_array[index % LOG_SIZE].counter, 1, MPI_INT, probe_status.MPI_SOURCE , DATA_TAG, slave_comm, &respond_request);
                    // wait to make sure send is done before testing for end cond, otherwise a process might be left hanging
                    if(exit_all)
                        MPI_Wait( &respond_request , &respond_status);
                }
                index++; 
            }
            // barrier to ensure that all threads are in sync before moving on to the next simulated time iteration
            #pragma omp barrier
        }
    }
    // finally send a summary report to the base station before terminating.
    struct Summary summary;
    summary.total_master_msgs = master_message_count;
    summary.total_nbr_msgs = neighbour_message_count;
    MPI_Send( &summary, 1, mpi_summary , base_station, SUMMARY_TAG , master_comm);
    //MPI_Barrier( slave_comm);
    MPI_Comm_free( &slave_comm );
    fclose(pFile);  
}