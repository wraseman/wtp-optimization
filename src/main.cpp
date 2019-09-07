/* main.cpp */

#include "wtp.h"
#include "auto_dose.h"
#include "wtp_optimize.h"

// Global variables
std::string current_datetime;
SimOptParameters simopt_params;

int main(int argc, char *argv[])
{
    /* Record current date and time */
    time_t rawtime;
    struct tm *info;
    char buffer[256];

    time(&rawtime);
    info = localtime(&rawtime);
    strftime(buffer, sizeof(buffer), "%m-%d-%y_%H-%M-%S", info);
    current_datetime = buffer; /* String which contains the date and time at the start of running the program */

    /* Read in water treatment plant model input file */
    ProcessTrain *train;
    if ((train = AllocProcessTrain()) == NULL) // Allocate memory for process train pointer for WTP Model
    {
        fprintf(stderr, "Error: cannot allocate memory for process train.\n");
    }

    const char* wtp_filepath= "./in/wtp_train/conv.wtp";
    open_wtp(wtp_filepath, train, NULL);             // Open and read in water treatment data
    // open_wtp(wtp_filepath, train, stdout); 

    // Initialize variables for influent and operations data
    std::string in_dir = "./in/single_sim/"; // input directory
    std::string influent_dir = "influent/";  // influent directory
    std::string operations_dir = "operations/";  // operations directory

    struct InfluentParameters *influent_wq;
    std::vector<double> influent_data;
    std::string influent_file;
    std::string influent_filepath;

    struct OperationalParameters *operations;
    std::vector<double> operations_data;
    std::string operations_file;  // operations file
    std::string operations_filepath;

    // Initialize variables for optimization run mode
    int num_func_evals; 
    int num_wq_scenarios; 

    // Read and validate command line arguments 
    int opt;
    std::string runmode;
    int f_flag = FALSE;
    int n_flag = FALSE;
    int i_flag = FALSE;
    int o_flag = FALSE;

    while ((opt = getopt(argc, argv, "r:f:n:i:o:h")) != -1)
    {
        switch (opt)
        {
        case 'r': // run mode: either simulate or optimize
            runmode = std::string(optarg);
            validate_optarg_runmode(runmode); // check that run mode is valid
            printf("run mode: %s\n", runmode.c_str());
            break;

        case 'f':                             // number of function evalutions (optimization mode only)
            validate_optarg_int(optarg, opt); // make sure input is non-zero integer
            printf("number of function evaluations: %s\n", optarg);
            num_func_evals = atoi(optarg); 
            f_flag = TRUE;
            break;

        case 'n':                             // number of influent scenarios (optimization mode only)
            validate_optarg_int(optarg, opt); // make sure input is non-zero integer
            printf("number of influent scenarios: %s\n", optarg);
            num_wq_scenarios = atoi(optarg); 
            n_flag = TRUE;
            break;

        case 'i':                              // influent file (simulation mode only)
            // validate_optarg_file(influent_filepath.c_str(), opt); // check that file exists
            printf("influent file: %s\n", optarg);
            i_flag = TRUE;

            // Read influent data
            influent_file = optarg;  // input file name
            influent_filepath = in_dir + influent_dir + influent_file;
            influent_data = read_single_sim(influent_filepath, true);  // read csv file with operations data

            // Populate influent parameter structure
            influent_wq = new InfluentParameters;
            influent_wq->alkalinity = influent_data[INF_ALKALINITY];
            influent_wq->nh3 = influent_data[INF_AMMONIA];
            influent_wq->bromide = influent_data[INF_BROMIDE];
            influent_wq->calcium = influent_data[INF_CALCIUM_HARDNESS];
            influent_wq->hardness = influent_data[INF_TOTAL_HARDNESS];
            influent_wq->pH = influent_data[INF_PH];
            influent_wq->temp = influent_data[INF_TEMPERATURE];
            influent_wq->toc = influent_data[INF_TOTAL_ORGANIC_CARBON];
            influent_wq->ntu = influent_data[INF_TURBIDITY];
            influent_wq->uv254 = influent_data[INF_UV254];
            break;

        case 'o':                              // operations file (simulation mode only)
            // validate_optarg_file(optarg, opt); // check that file exists
            printf("operations file: %s\n", optarg);
            o_flag = TRUE;

            // Read influent data
            operations_file = optarg;  // operations file
            operations_filepath = in_dir + operations_dir + operations_file;
            operations_data = read_single_sim(operations_filepath, true);  // read csv file with operations data

            operations = new OperationalParameters;
            operations->alk_setpt = operations_data[0];
            operations->pH_setpt = operations_data[1];
            operations->DBPsf = operations_data[2];
            break;

        case 'h': // help
            display_usage_help();
            break;

        default : break;
        
        }
    }

    // Run simulation or optimization
    if (runmode.compare("simulate") == 0)  // simulate mode
    {
        // Validate simulate parameters
        if ((i_flag == FALSE) || (o_flag == FALSE)) 
        {
            fprintf(stderr, "Error: both influent and operations files must be specified using the \"-i\" and \"-o\" flags.\nFor additional documentation, use \"-h\".\n");
            exit(EXIT_FAILURE);
        }

        /* Record command line arguments (i.e., run configuration) */
        std::string cli_dir; // command line interface directory
        std::string cli_ext; // command line interface extension
        cli_dir = "./out/single_sim/cli_args/";
        cli_ext = ".cli";
        std::string filepath_cli_args = cli_dir + current_datetime + cli_ext;

        save_cli_args(filepath_cli_args, argc, argv);

        /* Copy .wtp file used for the simulation */
        std::string copy_dir = "./out/single_sim/wtp_train/"; 
        std::string copy_ext = ".wtp";
        std::string copy_filepath = copy_dir + current_datetime + copy_ext;
        copy_file(wtp_filepath, copy_filepath);

        /* Define output file for simulation results */
        std::string result_dir; // result directory
        std::string result_ext; // result extension
        result_dir = "./out/single_sim/result/";
        result_ext = ".result";
        std::string filepath_result = result_dir + current_datetime + result_ext;

        FILE *fres = fopen(filepath_result.c_str(), "w"); // open file stream

        if (!fres)
        { /* Error handling */
            printf("Error opening %s \n, check that the proper directory and file name has been specified.\n", filepath_result.c_str());
            exit(EXIT_FAILURE);
        }

        std::cout << "Single simulation with automatic chemical dosing" << std::endl;

        /* Call single simulation with automatic dosing */
        // single_sim(train, current_datetime, fres);
        single_sim_mode(train, operations, influent_wq, fres); // TODO: call single_sim() instead for continuity

        fclose(fres); // close file stream

        // Free dynamically allocated memory
        delete operations;
        delete influent_wq;
    }
    else if (runmode.compare("optimize") == 0)  // optimize mode
    {
        // Validate optimize parameters
        if ((f_flag == FALSE) || (n_flag == FALSE))
        {
            fprintf(stderr, "Error: both number of function evaluations and influent scenarios must be specified using the \"-f\" and \"-n\" flags.\nFor additional documentation, use \"-h\".\n");
            exit(EXIT_FAILURE);
        }

        /* Record command line arguments */
        std::string cli_dir; // command line interface directory
        std::string cli_ext; // command line interface extension
        cli_dir = "./out/sim_opt/cli_args/";
        cli_ext = ".cli";
        std::string filepath_cli_args = cli_dir + current_datetime + cli_ext;

        save_cli_args(filepath_cli_args, argc, argv);

        /* Copy .wtp file used for the optimization */
        std::string copy_dir = "./out/sim_opt/wtp_train/"; 
        std::string copy_ext = ".wtp";
        std::string copy_filepath = copy_dir + current_datetime + copy_ext;
        copy_file(wtp_filepath, copy_filepath);

        // Optimize water treatment
        simopt_params.num_func_evals = num_func_evals;
        simopt_params.num_wq_scenarios = num_wq_scenarios;

        // /* Define Monte Carlo parameters for influent water quality data */
        // params.mc.mc_flag = true;
        // params.mc.n_samples = num_wq_scenarios;

        // /* Define time series parameters for influent water quality data */
        // params.ts.ts_flag = true;
        // params.ts.n_years = 11;           
        // params.ts.timesteps_per_year = 4; // quarterly data

        // /* Define optimization problem */
        // params.num_func_evals = num_func_evals;
        // params.problem = 1;

        std::cout << "Simulation-optimization mode\n";
        sim_opt_mode(train);
    }

    /* Free memory */
    FreeProcessTrain(train);

    /* Pause to read output */
    // system("pause");

    return 0;
}

/* Validation functions */
// purpose: validate run mode command line argument is either "simulate" or "optimize" (optimize)
void validate_optarg_runmode(std::string runmode)
{
    if ((runmode.compare("simulate") == 0) || (runmode.compare("optimize") == 0))
    {
        // If true, then entered run mode is either "simulate" or "optimize". Do nothing.
    }
    else
    { // If false, an invalid run mode has been entered.
        fprintf(stderr, "Error: \"simulate\" and \"optimize\" are the only valid run mode options\n");
        exit(EXIT_FAILURE);
    }
    return;
}

// purpose: validate that command line argument is a non-zero integer
void validate_optarg_int(char *optarg, char opt)
{
    if (atoi(optarg) == 0)
    {
        fprintf(stderr, "Error: value for \"-%c\" argument must be a non-zero integer\n", opt);
        exit(EXIT_FAILURE);
    }
    return;
}

// purpose: check that file specified in command line argument exists
void validate_optarg_file(const char *optarg, char opt)
{
    FILE *file;
    if ((file = fopen(optarg, "r")))
    {
        // file exists
        fclose(file);
    }
    else
    {
        // file does not exist
        fprintf(stderr, "Error: the file \"%s\" specified in \"-%c\" argument does not exist\n", optarg, opt);
        exit(EXIT_FAILURE);
    }
    return;
}

/* Documentation function */
// purpose: display command line arguments documentation
void display_usage_help()
{
    printf("\nCommand line arguments available for wtp-optimize:\n");
    printf("-r (run mode): enter either \"simulate\" or \"optimize\"\n");
    printf("-f (function evalutions): enter a non-zero integer [optimization mode only]\n");
    printf("-n (number of influent scenarios): enter a non-zero integer [optimization mode only]\n");
    printf("-i (influent file): enter relative path to influent directory [simulate mode only]\n");
    printf("-o (operations file): enter relative path to operations directory [simulate mode only]\n");
    printf("-h (help): display command line arguments documentation\n");
    printf("\n");
    printf("Simulation example (if executable is in the binary directory):\n");
    printf("./bin/wtp-optimize.exe -r simulate -i ./in/influent.txt -o ./in/operations.txt\n");
    printf("\n");
    printf("Optimization example (if executable is in the binary directory):\n");
    printf("./bin/wtp-optimize.exe -r optimize -f 10000 -n 100\n");
    return;
}

/* Reproducibility functions */
// purpose: save file of command line arguments used for each simulation or optimization run
void save_cli_args(std::string filepath_cli_args, int argc, char** argv)
{
        FILE *fcon = fopen(filepath_cli_args.c_str(),"a");

        if (!fcon)
        { /* Error handling */
            printf("Error opening %s \n, check that the proper directory and file name has been specified.\n", filepath_cli_args.c_str());
            exit(EXIT_FAILURE);
        }
        fprintf(fcon, "Command line arguments used to configure simulation:\n");
        for (int i = 0; i < argc; ++i) 
        {
            fprintf(fcon, "%s ", argv[i]);
        }

        fclose(fcon);

        return;
}

// purpose: copy a file from a source location to another destination
// source: https://www.geeksforgeeks.org/c-program-copy-contents-one-file-another-file/ (edited)
void copy_file(std::string source_path, std::string dest_path) 
{ 
    FILE *fptr1, *fptr2; 
    char c; 
  
    // Open one file for reading 
    fptr1 = fopen(source_path.c_str(), "r"); 
    if (fptr1 == NULL) 
    { 
        printf("Cannot open file %s \n", source_path.c_str()); 
        exit(0); 
    } 
  
    // Open another file for writing 
    fptr2 = fopen(dest_path.c_str(), "w"); 
    if (fptr2 == NULL) 
    { 
        printf("Cannot open file %s \n", dest_path.c_str()); 
        exit(0); 
    } 
  
    // Read contents from file 
    c = fgetc(fptr1); 
    while (c != EOF) 
    { 
        fputc(c, fptr2); 
        c = fgetc(fptr1); 
    } 
  
    fclose(fptr1); 
    fclose(fptr2); 

    return;
}

// int main()
// {
//     /* Record current date and time */
//     time_t rawtime;
//     struct tm *info;
//     char buffer[256];

//     time(&rawtime);
//     info = localtime(&rawtime);
//     strftime(buffer, sizeof(buffer), "%m-%d-%y_%H-%M-%S", info);
//     current_datetime = buffer; /* String which contains the date and time at the start of running the program */

//     // TODO: read this in from JSON instead
//     int run_mode = 2;
//     // int run_mode = 4;  // simulation optimizatino run

//     // /* Read JSON run configuration file */
//     // const char *config_filename = "config.json";
//     // std::string str = read_json(config_filename);
//     // const char *json = str.c_str(); // convert string to c string (const char *)

//     // /* Parse JSON data */
//     // rapidjson::Document document;
//     // document.Parse(json);

//     /* Read in treatment type data from configuration file and create associated WTP treatment train */
//     ProcessTrain *train;
//     if ((train = AllocProcessTrain()) == NULL) // Allocate memory for process train pointer for WTP Model
//     {
//         fprintf(stderr, "Error: cannot allocate memory for process train.\n");
//     }

//     // TODO: read this in from JSON instead
//     const char* wtp_filename = "./in/wtp_train/conv.wtp";
//     open_wtp(wtp_filename, train, stdout);             // Open and read in water treatment data

//     // const char *treat_name = "treatment type";
//     // if (document.HasMember(treat_name) && document[treat_name].IsString())
//     // {
//     //     std::string treat_type = document[treat_name].GetString(); // Get string of treatment type
//     //     std::string wtp_filename = get_wtp_filename(treat_type);   // Get treatment train information from WTP Model input file
//     //     open_wtp(wtp_filename.c_str(), train, stdout);             // Open and read in water treatment data
//     // }
//     // else
//     // {
//     //     fprintf(stderr, "Error: treatment type either not specified or not of type string in config.json");
//     //     exit(EXIT_FAILURE);
//     // }

//     // /* Read in run mode data from configuration file */
//     // int run_mode;
//     // const char *mode_name = "run mode";
//     // if (document.HasMember(mode_name) && document[mode_name].IsInt())
//     // {
//     //     run_mode = document[mode_name].GetInt();
//     // }
//     // else
//     // {
//     //     fprintf(stderr, "Error: run mode either not specified or not of type integer in config.json");
//     //     exit(EXIT_FAILURE);
//     // }

//     /* Copy configuration file for reference */
//     // source: https://stackoverflow.com/questions/10195343/copy-a-file-in-a-sane-safe-and-efficient-way
//     std::string config_dir;
//     if (run_mode == 2)
//     {  // single simulation
//         config_dir = "./out/single_sim/config/";
//     }
//     else if (run_mode == 4)
//     {  // simulation-optimization
//         config_dir = "./out/sim_opt/config/";
//     }
//     else
//     {
//         fprintf(stderr, "Error: configuration file not supported for this run mode.");
//         exit(EXIT_FAILURE);
//     }
//     std::string config_name = config_dir + current_datetime + ".json";
//     std::ifstream src("config.json", std::ios::binary);
//     std::ofstream dst(config_name.c_str(), std::ios::binary);  // need to convert inputfile to c string for code that is C++03 or below (source: https://stackoverflow.com/questions/35045781/no-matching-function-for-call-to-stdbasic-ofstreamchar-stdchar-traitscha)
//     dst << src.rdbuf();  // copy configuration file to destination folder

//     /* Based on run mode, call appropriate simulation/optimization function */
//     switch (run_mode)
//     {
//     // TODO: remove this case. Don't need to run simple WTP Model run using this program
//     // case 1:
//     // {
//     //     std::cout << "Single simulation mode\n";

//     //     // Read in single simulation run parameters
//     //     // TODO: read in parameters (don't just hard code "true" in populate run parameter structure)
//     //     /* Read in run mode data from configuration file */
//     //     // const char *treat_name = "treatment type";
//     //     // bool treat_hasmember = document.HasMember(treat_name);
//     //     // bool treat_isstring = document[treat_name].IsString();
//     //     // string treat_getstring = document[treat_name].GetString();
//     //     // string treatment_type = treat_type_config(treat_hasmember, treat_isstring, treat_getstring);

//     //     // Populate run parameter structure
//     //     SingleSimParameters params;
//     //     params.list_process_train = true;
//     //     params.run_model = true;
//     //     params.thm_and_disinfection = true;

//     //     single_sim_mode(train, params);

//     //     /* Copy configuration file for reference */
//     //     // TODO: copy configuration file!
//     //     // std::filesystem::copy("config.json", "./out/single_sim")

//     //     break;
//     // }
//     case 2:
//     {
//         // Read in single simulation run parameters
//         // TODO: are there parameters that need to be read in?
//         // TODO: read in parameters

//         /* Read in run parameters from JSON */
//         // TODO: read these in from config.json

//         /* Read in single simulation run operations and influent water quality */
//         std::vector<double> operations_data;
//         std::vector<double> influent_data;
//         std::string in_dir = "./in/single_sim/"; // input directory

//         // std::string operations_file = "/operations/2019-03-04_operations.csv";  // operations file
//         std::string operations_file = "/operations/CLP-Q4_cluster-3_operations.csv";  // operations file
//         std::string operations_filepath = in_dir + operations_file;
//         operations_data = read_single_sim(operations_filepath, true);  // read csv file with operations data

//         struct OperationalParameters *operations;
//         operations = new OperationalParameters;
//         operations->alk_setpt = operations_data[0];
//         operations->pH_setpt = operations_data[1];
//         operations->DBPsf = operations_data[2];

//         std::string influent_file = "/influent/CLP-9999-Q4_influent.csv";  // input file name
//         std::string influent_filepath = in_dir + influent_file;
//         influent_data = read_single_sim(influent_filepath, true);  // read csv file with operations data

//         struct InfluentParameters *influent_wq;
//         influent_wq = new InfluentParameters;
//         influent_wq->alkalinity = influent_data[INF_ALKALINITY];
//         influent_wq->nh3 = influent_data[INF_AMMONIA];
//         influent_wq->bromide = influent_data[INF_BROMIDE];
//         influent_wq->calcium = influent_data[INF_CALCIUM_HARDNESS];
//         influent_wq->hardness = influent_data[INF_TOTAL_HARDNESS];
//         influent_wq->pH = influent_data[INF_PH];
//         influent_wq->temp = influent_data[INF_TEMPERATURE];
//         influent_wq->toc = influent_data[INF_TOTAL_ORGANIC_CARBON];
//         influent_wq->ntu = influent_data[INF_TURBIDITY];
//         influent_wq->uv254 = influent_data[INF_UV254];

//         /* Define output file for simulation results */
//         std::string result_dir; // result directory
//         std::string result_ext; // result extension
//         result_dir = "./out/single_sim/result/";
//         result_ext = ".result";
//         std::string filename_result = result_dir + current_datetime + result_ext;

//         FILE *fres = fopen(filename_result.c_str(), "w"); // open file stream

//         if (!fres)
//         { /* Error handling */
//             printf("Error opening %s \n, check that the proper directory and file name has been specified.\n", filename_result.c_str());
//             exit(EXIT_FAILURE);
//         }

//         std::cout << "Single simulation with automatic chemical dosing" << std::endl;

//         /* Call single simulation with automatic dosing */
//         // single_sim(train, current_datetime, fres);
//         single_sim_mode(train, operations, influent_wq, fres); // TODO: call single_sim() instead for continuity

//         fclose(fres); // close file stream

//         // Free dynamically allocated memory
//         delete operations;
//         delete influent_wq;

//         break;
//     }
//     // case 3:
//     // {
//     //     cout << "Multi-simulation mode";
//     //     multi_sim_mode(train, params);
//     //     break;
//     // }
//     case 4:
//     {
//         // Read in simulation-optimization run parameters
//         // TODO: read in simulation-optimization parameters

//         SimOptParameters params;

//         /* Define Monte Carlo parameters for influent water quality data */
//         params.mc.mc_flag = true;
//         params.mc.n_samples = 5;

//         /* Define time series parameters for influent water quality data */
//         params.ts.ts_flag = true;
//         params.ts.n_years = 5;            // TODO: check if this the correct number of years
//         params.ts.timesteps_per_year = 4; // quarterly data

//         /* Define optimization problem */
//         params.problem = 1;

//         std::cout << "Simulation-optimization mode\n";
//         sim_opt_mode(train, params);
//         break;
//     }
//     // case 5:
//     // {
//     //     cout << "Validation mode";
//     //     validation_mode(train, params);
//     //     break;
//     // }
//     default:
//     {
//         fprintf(stderr, "Error: run mode not supported. Please specify another run mode.");
//         exit(EXIT_FAILURE);
//         break;
//     }
//     }

//     /* Free memory */
//     FreeProcessTrain(train);

//     /* Pause to read output */
//     system("pause");

//     return 0;
// }

// /* Declare main functions */

// /* Read contents of run configuration file (config.json) into a string */
// std::string read_json(const char *filename)
// {
//     // source: https://stackoverflow.com/questions/2602013/read-whole-ascii-file-into-c-stdstring
//     // std::ifstream file(filename);
//     // file.seekg(0, std::ios::end);
//     // size_t size = file.tellg();
//     // std::string buffer(size, ' ');
//     // file.seekg(0);
//     // file.read(&buffer[0], size); // read contents into buffer
//     // file.close();                // close input file stream

//     std::ifstream t(filename);
//     std::stringstream buffer;
//     buffer << t.rdbuf();

//     return buffer.str();
// }

// /* Return WTP Model input file name based on config.json treatment type */
// std::string get_wtp_filename(std::string treatment_type)
// {
//     std::string wtp_filename;

//     if (treatment_type == "conventional")
//     {
//         wtp_filename = "./in/wtp_train/conv.wtp";
//     }
//     else // Unsupported treatment type
//     {
//         fprintf(stderr, "Error: only conventional treatment plant supported. config.json treatment type must be \"conventional\".");
//         exit(EXIT_FAILURE);
//     }

//     return wtp_filename;
// }