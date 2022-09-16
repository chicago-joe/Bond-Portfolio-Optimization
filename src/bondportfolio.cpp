

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cstdlib>

#include "lp_lib.h"

using namespace std;

const double eRROR = 1e-10; // for use in Newton Raphson Method - this is a very small number that is > zero
int bonds;  // renamed for # of bonds rather than cashflows. in Finance jargon: each bond has its own # of cashflows
vector <double> current_price_list;     // current prices of bonds-> note: name is interchangeable with PV of Bond
vector <int> maturity_list;             // stores bond maturities (ie # of cash flows of each bond)
vector <double> yield_to_maturity;
vector <double> duration;
vector <double> convexity;
double investment_amt;                  // renamed for portfolios rather than paying off a debt obligation
double investment_horizon;              // renamed for portfolios, this is "time that debt obligation is due"
vector <double> req_pct_of_bond;        // following the above, this is the required % of the FV of each bond

vector <vector <double>> cash_flow_list;

double function(vector <double> cash_flow, double price, int maturity, double rate)
{
    // computes the ytm of bond that will equate the bond's value with the sum of the FV of its cashflow pmts
    // hence, it is named discount r for "discount rate"
    double discountr = price * pow((1 + rate), maturity);
    for (int i = 0; i < maturity; i++)
    {
        discountr -= (cash_flow[i] * pow((1 + rate), (maturity - 1 - i)));
    }
    return (discountr);
}

double derivative_function(vector <double> cash_flow, double price, int maturity, double rate)
{
    // computes f'(r) for use in the Newton_Raphson function of the next section
    double fprime = maturity * price * pow((1 + rate), maturity - 1);
    for (int i = 0; i < maturity - 1; i++)
    {
        fprime -= ((maturity - 1 - i) * cash_flow[i] * pow((1 + rate), (maturity - 2 - i)));
    }
    return (fprime);
}

double Newton_Raphson(vector <double> cash_flow, double price, int maturity, double rate)
{
    // write a function that finds the (only) +ve root of f(r) of page 2 of lesson 3 using Newton-Raphson method
    double discountr, fprime, root_r = 0.5;

    discountr = function(cash_flow, price, maturity, root_r);
    fprime = derivative_function(cash_flow, price, maturity, root_r);

    // as mentioned in line 14, this is an extremely small number close to 0.
    // this solves the flaw that occurs when iterating NewtownRaphson
    // ie NewtonRaphson doesn't fail with strictly +cashflows
    while (abs(discountr) > eRROR)
    {
        discountr = function(cash_flow, price, maturity, root_r);
        fprime = derivative_function(cash_flow, price, maturity, root_r);
        root_r -= discountr / fprime;
    }
    return (root_r);
}

double get_duration(vector <double> cash_flow, double price, int maturity, double rate)
{
    // computes each bond's duration
    // duration = first derivative of bond price wrt (discount rate / bond price)
    double duration = 0.0;
    double delta = 0.0;

    for (int i = 0; i < maturity; i++)
    {
        // Vb + change_in_Vb ~ delta --> approximately. (this note can be ignored, purpose was for pseudo-code)
        delta += (cash_flow[i] * (i + 1) / pow((1 + rate), i + 1));
    }
    duration = delta / price;
    return (duration);
}

double get_convexity(vector <double> cash_flow, double price, int maturity, double rate)
{
    // computes each bond's convexity
    // convexity = second-derivative of bond price wrt (discount rate / bond price)
    double convexity = 0.0;
    double delta = 0.0;

    for (int i = 0; i < maturity; i++)
    {
        delta += (cash_flow[i] * ((i + 1)*(i + 2)) / pow((1 + rate), i + 3));
    }
    convexity = delta / price;
    return(convexity);
}

double present_value_of_debt()
{
    // compute PV of future debt obligation using the average-value-of-the-YTMs
    double avg_ytm = 0.0;
    double sum_ytm = 0.0;

    for (int i = 0; i < bonds; i++)
    {
        // sum the YTMs
        sum_ytm += yield_to_maturity[i];
    }
    // compute the average YTM = sumYTM / #ofbonds
    avg_ytm = sum_ytm / bonds;
    return (investment_amt / (pow((1 + avg_ytm), investment_horizon)));
}

// Replica of part of present_value_of_debt calculation
// to be called when outputting "Average YTM (to compute PV of debt)"
double average_ytm()
{
    double sum_ytm = 0.0;
    for (int i = 0; i < bonds; i++)
    {
        sum_ytm += yield_to_maturity[i];
    };
    return (sum_ytm / bonds);
}

void print_data(char *filename)
{
    cout << "Input File: " << filename << endl << endl;
    cout << "We want to invest " << investment_amt << " over " << investment_horizon << " years" << endl;
    cout << "Number of Bonds: " << bonds << endl;
    for (int i = 0; i < bonds; i++)
    {
        cout << "------------------------------------------------------" << endl;
        cout << "Bond #" << i + 1 << endl;
        cout << "Current Price = " << current_price_list[i] << endl;
        cout << "Maturity = " << maturity_list[i] << endl;
        cout << "Percentage of Face Value that would meet investment requirements = " << req_pct_of_bond[i] << endl;
        cout << "Yield to Maturity = " << yield_to_maturity[i] << endl;

        // LP-model formulation takes into account the required % of each bond needed.
        // the following statements show how we incorporate the % requirements of each bond
        // when using lpsolve to meet constraints and find the optimal portfolio

        cout << "Duration = " << duration[i] << endl;   // (prior to dividing by % of face value)
        cout << "Duration (to be used in LP-formulation below) = " << duration[i] / req_pct_of_bond[i] << endl;
        cout << "(Note) " << duration[i] << " = " << duration[i] / req_pct_of_bond[i] << " x "
            << req_pct_of_bond[i] << endl;
        cout << "Convexity = " << convexity[i] << endl;     // (prior to dividing by % of face value)
        cout << "Convexity (to be used in LP-formulation below) = " << convexity[i] / req_pct_of_bond[i] << endl;
        cout << "(Note) " << convexity[i] << " = " << convexity[i] / req_pct_of_bond[i] << " x "
            << req_pct_of_bond[i] << endl;
    }
    cout << "------------------------------------------------------" << endl;
    cout << "Average YTM (to compute PV of debt) = " << average_ytm() << endl;
    cout << "Present value of debt = " << present_value_of_debt() << endl;
    cout << "------------------------------------------------------" << endl;
}

void get_data(char* argv[])
{
    // write the code that reads the data from the file identified on the command-line.
    ifstream input_file(argv[1]);

    // verify input file is in local directory
    if (input_file.is_open())
    {
        input_file >> bonds;
        for (int i = 0; i < bonds; i++)
        {
            vector <double> cash_flow;
            double cf; 						    // temp for cashflows
            double p; 						    // temp for current price
            int m;							    // temp for maturity (ie investment horizon)

            input_file >> p;				    // read price from file
            current_price_list.push_back(p);    // push price into current_price_list

            input_file >> m;				    // read maturity from input file
            maturity_list.push_back(m);		    // push maturity into maturity_list

            for (int k = 0; k <= m - 1; k++)
            {
                input_file >> cf;               // same as above -> read file and push into list
                cash_flow.push_back(cf);        // same as above -> read file and push into list
            }
            // push CFs into cash_flow_list
            cash_flow_list.push_back(cash_flow);

            // compute yield to maturity of bond (ytm) using NewtonRaphson Method
            double r = Newton_Raphson(cash_flow, p, m, 0.0);
            yield_to_maturity.push_back(r);


            double d = get_duration(cash_flow, p, m, r);  // get duration of bond
            duration.push_back(d);                        // push into list

            // get convexity of the cash flow
            double c = get_convexity(cash_flow, p, m, r);  // get duration of bond
            convexity.push_back(c);                        // push into list
        }
        // reads investment amount (obligation) and time horizon (time obligation is due) from file
        input_file >> investment_amt;
        input_file >> investment_horizon;

        for (size_t i = 0; i < bonds; i++)
        {
            // computing the percentage of each bond needed to satisfy investment requirements
            req_pct_of_bond.push_back(present_value_of_debt() / current_price_list[i]);
        }
    }
    else
    {
        // if input file not found in local directory...
        cout << "ERROR: FILE DOES NOT EXIST IN LOCAL DIRECTORY." << endl;
        system("pause");
    }
}

void get_optimal_portfolio()
{
    // write the lp_solve specific material that computes the optimal_portfolio
    lprec *lp;
    // make a 1d array
    double* solution = new double[bonds];
    double* row = new double[bonds + 1];
    row[0] = 0;

    lp = make_lp(0, bonds);
    set_verbose(lp, 3);

    for (int i = 0; i < bonds; i++)
    {
        row[i + 1] = current_price_list[i];
    }

    // Adds constraint R1 that solution amount (value) must = the present value of debt
    add_constraint(lp, row, EQ, present_value_of_debt());

    for (int i = 0; i < bonds; i++)
    {
        row[i + 1] = (duration[i] / req_pct_of_bond[i]);
    }
    // Adds constraint R2 that solution duration must = the investment horizon (maturity)
    add_constraint(lp, row, EQ, investment_horizon);

    for (int i = 0; i < bonds; i++)
    {
        row[i + 1] = convexity[i] / (req_pct_of_bond[i]);
    }

    set_maxim(lp); // used the maximize function of lpsolve --> same as -(minimize lp)
    set_obj_fn(lp, row); // sets convexity objective of lpsolve --> objective: maximize convexity
    print_lp(lp);

    if (solve(lp))
    {
        // if there is no portfolio that will meet our investment requirements (aka obligation)
        // and immunize against parallel shifts in the term structure...
        cout << "There is no portfolio that meets the duration constraint of " << investment_horizon << " years" << endl;
    }
    else
    {
        get_variables(lp, solution);
        // shows lp objective: portfolio w/ largest convexity among all possible portfolio choices
        cout << "Largest convexity we can get is: " << get_objective(lp) << endl << endl << "Optimal portfolio (%): " << endl;
        for (int i = 0; i < bonds; i++)
        {
            cout << "Bond " << i + 1 << ": " << solution[i] << endl;
        }
     
        cout << endl << "To create optimal portfolio: BUY " << endl;
        for (int i = 0; i < bonds; i++)
            if (solution[i] > 0)
                // solution: optimal portfolio is to buy: lambda(i)% of bond(i)..
                cout << "$" << solution[i] * current_price_list[i] << " of Bond " << i + 1 << endl;

    }


    delete_lp(lp);

    delete row;
    delete solution;
}

int main(int argc, char* argv[])
{
    if (argc == 1) {
        cout << "ERROR: INPUT FILE NOT FOUND. " << endl;
    }
    else
    {
        get_data(argv);

        print_data(argv[1]);

        get_optimal_portfolio();
        cout << endl << endl;
        system("pause");
    }
    return (0);
}

