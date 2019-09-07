#include <string>
#include <vector>
#include <sstream>  //istringstream
#include <iostream> // cout
#include <fstream>  // ifstream
#include <cstdlib>
#include <stdexcept>      // std::invalid_argument

using namespace std;

/* read_csv.cpp */
/* Purpose: this file contains functions which support reading .csv (comma separated value) input
    files for different run modes. This includes single and multi-simulation modes and Monte Carlo
    influent scenarios. 
*/

// define macros about Monte Carlo simulation of influent water quality input data
#define N_YEARS 11            // number of years on record
#define N_QUARTERS_PER_YEAR 4 // number of quarters in a year
#define MAX_N_SIMS 5000       // maximum number of simulations that the user can run

/**
 * Reads csv file for single simulation mode containing decision variables for a single quarter. 
 * source: https://waterprogramming.wordpress.com/2017/08/20/reading-csv-files-in-c/
 * @param inputFileName: input file name (full path).
 * @param header: skip 1st row of text which contains column names? (true/false)
 * @return data as vector of doubles.
 */

vector<double> read_single_sim(string filename, bool header)
{
    std::ifstream file(filename.c_str());
    vector<double> data;
    string str;  // string read from file

    /* Define maximum number of lines that can be read */
    int max_nlines = 1;
    if (header == true)
    {
        max_nlines++;
    }

    /* Error handling */
    if (!file)
    {
        cerr << "Could not read file " << filename << "\n";
        exit(EXIT_FAILURE);
    }

    /* Define number of lines to be read based on number of simulations */
    int nlines = max_nlines;
    int line = 0; // initialize current line number
    string buffer;

    while (line <= nlines)
    {
        line++; // add count to current line number

        // Read each line of text into a buffer
        if (!getline(file, buffer))
        {
            break;
        }
        istringstream stream_buffer(buffer);

        // If there is a header, skip the first line
        if (header == true && line == 1)
        {
            continue;
        }

        // Read each comma-separated value of the line
        while (stream_buffer)
        {
            string value_str;
            if (!getline(stream_buffer, value_str, ','))
            {
                break;
            }
            try
            {
                data.push_back(stof(value_str));
            }
            catch (const std::invalid_argument e)
            {
                cout << "NaN found in file " << filename << " line " << line << endl;
                e.what();
            }
        }
    }

    return data;
}

/**
 * Reads csv file into table, exported as a vector of vector of doubles.
 * source: https://waterprogramming.wordpress.com/2017/08/20/reading-csv-files-in-c/
 * @param inputFileName: input file name (full path).
 * @param header: skip 1st row of text which contains column names? (true/false)
 * @return data as vector of vector of doubles.
 */

vector< vector<double> > read_montecarlo(string filename, bool header, int nsims)
{

    vector< vector<double> > data;
    std::ifstream file(filename.c_str());
    
    string str;  // string read from file

    /* Define maximum number of lines that can be read */
    int max_nlines = N_YEARS * N_QUARTERS_PER_YEAR * MAX_N_SIMS;
    if (header == true)
    {
        max_nlines++;
    }
    /* Error handling */
    if (nsims > MAX_N_SIMS)
    {
        cerr << "Number of Monte Carlo simulations specified (" << nsims
             << ") exceeds maximum allowable number (" << MAX_N_SIMS << ")"
             << "\n";
        exit(EXIT_FAILURE);
    }

    /* Error handling */
    if (!file)
    {
        cerr << "Could not read file " << filename << "\n";
        exit(EXIT_FAILURE);
    }

    /* Define number of lines to be read based on number of simulations */
    int nlines = N_YEARS * N_QUARTERS_PER_YEAR * nsims;
    int line = 0; // initialize current line number
    string buffer;

    while (line <= nlines)
    {
        line++; // add count to current line number

        // Read each line of text into a buffer
        if (!getline(file, buffer))
        {
            break;
        }
        istringstream stream_buffer(buffer);
        vector<double> record;

        // If there is a header, skip the first line
        if (header == true && line == 1)
        {
            continue;
        }

        // Read each comma-separated value of the line
        while (stream_buffer)
        {
            string value_str;
            if (!getline(stream_buffer, value_str, ','))
            {
                break;
            }
            try
            {
                record.push_back(stof(value_str));
            }
            catch (const std::invalid_argument &e)
            {
                cout << "NaN found in file " << filename << " line " << line << endl;
                e.what();
            }
        }

        data.push_back(record);
    }

    return data;
}