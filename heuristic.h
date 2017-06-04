#ifndef COIOTE_HEURISTIC_HEURISTIC_H
#define COIOTE_HEURISTIC_HEURISTIC_H

#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <iostream>
#include <list>
#include <map>
#include <thread>

using namespace std;

struct Data {
    /**
     * the costs ijkt
     */
    double**** costs;

    int ****app;

    /**
     * number of activity done by one user of type k
     */
    int* vetTaskOfPeople;  // numero di attività per ogni persona di un certo tipo

    /**
     * Activities to be done in node i during the given time horizon
     */
    int* vetTaskCells;  //numero di attività massime(da svolgere) in una cella

    /**
     * Number of users of type m in i during time period t (ikt)
     */
    int*** peopleTakenInCells;

    int **vetUsers; // contiene in ordine l'utente che costa meno


    /**
     * Vettore per ordinare le persone per division value
     */
    double ***vetDivision;

    int *objFunOfThreads;

    int **bestTypeTakenByThreads;



};

struct CellDestData {

    /**
     *  Division value per ogni cella
     */
    double* vetMinDivisione;
    /**
     *
     */
    int* vetUser;
    /**
     *
     */
    int* vetPeriod;
    /**
     *
     */
    int* vetCell;
    /**
     * Flag binario che mi dice se ho risolto una cella
     */
    int flag;
    /**
     * Dimensione del vettore di utenti ordinati per ogni cella
     */
    int size; //Dovrebbe partire logicamente da 1 ma a differenza del size classico lo facciam partire da 0 perché lo usiamo anche come indice
    /**
     * Numero della cella destinazione (per permettere il random)
     */
    int cellDest;
    /**
     * Numero di task da fare nella cella
     */
    int numberOfTask;

    /**
     * Task già eseguiti nella cella
     */
    int currentTask;
    /**
     * Indice da cui partire nel powerset
     */
    int startingPoint;
};

enum eFeasibleState {
    FEASIBLE,
    NOT_FEASIBLE_DEMAND,
    NOT_FEASIBLE_USERS
};

class Heuristic{
private:
    /**
     * Number of periods
     */
    int nTimeSteps;

    /**
     * Number of customer types
     */
    int nCustomerTypes;

    /**
     * Number of cells
     */
    int nCells;

    /**
     * Problem structure for parameters
     */
    Data problem;


    /**
     * Flag equals to true if the problems has a solution
     */
    bool hasSolution;

    /**
     * Variables of the problem (X in the model)
     */
    int**** solution;



public:
    /**
     * Default constructor
     * @return Heuristic object
     */
    Heuristic(){};

    /**
     * Constructor from external file
     * @param path path of the external file cotaining the instance of the problem
     * @return
     */
    Heuristic(string path);

    /**
     * Function to CHANGE!!! This function only makes a very bad solution for the problem
     * @param stat Array of statistics. In position 0 put the objVal, in position 1 the computational time
     * @param timeLimit Time limit for computation
     * @param verbose
     */
    void solveFast(vector<double>& stat, int timeLimit = - 1);

    /**
     * Puts KPIs in the statistics' array. Call this only if problem has a solution
     * @param stat Array of statistics
     */
    void getStatSolution(vector<double>& stat);

    /**
     * Write KPIs on a file
     * @param path path of the file
     * @param nameInstance name of the instance
     * @param stat array of statistics
     */
    void writeKPI(string path, string nameInstance, vector<double> stat);

    /**
     * Write the detailed solution on a file
     * @param path path of the file
     */
    void writeSolution(string path);

    /**
     * Check the feasibility of the problem
     * @param path path of the solution file
     * @return a state of the check (i.e. FEASIBLE if the solution is feasible)
     */

    eFeasibleState isFeasible(string path);

    void quickSort(double *arr, int *arr2 , int *arr3 , int *arr4 , int left, int right);
    void quickSort2(int *arr, int *arr2, int left,int right);

    void powerset(int *numberOfTask,int currentTask,int* vetUser,int* vetPeriod,
                  int* vetCell,int generalIndex,int *bestSolByType,
                  int* typeTaken,int *max,int curr,int cellDest,int ***pplTakenInCellCopy);

    void cellByCell(double tStart,int dimDestArray,int tId,int predefinedTaskCount,int cicloN);

    void democratic(double tStart,int dimDestArray,int tId,int predefinedTaskCount,int cicloN);


};

#endif //COIOTE_HEURISTIC_HEURISTIC_H
