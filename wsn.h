#define SHIFT_ROW 0
#define SHIFT_COL 1
#define DISP 1
#define MAX_ITERATION 50
#define NUM_PORTS 4
#define LOG_SIZE 8
#define EXIT_TAG 1
#define ASK_TAG 1000
#define DATA_TAG 2000
#define SUMMARY_TAG 3000
#define MAX_REPORTS 100
// Define a struct for datetime

struct DateTime {
    int year;
    int month;
    int day;
    int hour;
    int minute;
    int second;
};
struct Report{
    struct DateTime datetime;   // simulated datetime of the alert
    int messages_counter;       // counter for the number of messages sent to base station
    int nbrs_count;             // number of neighbours relative to the reporting node
    int reporting_node;         // rank of the reporting node
    int iteration;              // iteration in which the error occured
    double time_taken;          // time taken for communicating with neighbours
    double start_time;          // time at which the alert has been sent
    int nbrs[4]; // immediate neighbours: LEFT,RIGHT,UP,DOWN
};  

struct Summary{
    int total_nbr_msgs;         // total messages sent to neighbours before an alert
    int total_master_msgs;      // total messages sent to base station/master
};
// Define a struct for a log entry
struct Cycle {
    unsigned int cycle_num;     // the sumulated cycle number
    struct DateTime datetime;   // simulated datetime of the cycle (+10s each 0.5s)
    unsigned int counter;       // counter for the number of ports in use (0 means all ports are free to use)
};

