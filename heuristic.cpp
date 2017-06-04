#include <iostream>
#include <random>
#include "heuristic.h"
#include <thread>
#include <algorithm>
#include <cfloat>
#include "mingw.thread.h"

using namespace std;

Heuristic::Heuristic(string path){
    this->hasSolution = false;
    string line;
    string word;

    vector<float> vetProvaOrder;

    ifstream iffN(path.c_str());

    if (!iffN.is_open()) {
        cout << "Impossible to open" << path << endl;
        cin.get();
        exit(1);
    }

    getline(iffN, line);
    std::replace(line.begin(), line.end(), ';', ' ');
    istringstream iss(line);
    iss >> word;
    this->nCells = atoi(word.c_str());
    iss >> word;
    this->nTimeSteps = atoi(word.c_str());
    iss >> word;
    this->nCustomerTypes = atoi(word.c_str());

    // Memory allocation
    solution = new int***[nCells];
    problem.costs = new double***[nCells];
    for (int i = 0; i < this->nCells; i++) {
        problem.costs[i] = new double**[nCells];
        solution[i] = new int**[nCells];
        for (int j = 0; j < this->nCells; j++) {
            problem.costs[i][j] = new double*[nCustomerTypes];
            solution[i][j] = new int*[nCustomerTypes];
            for (int m = 0; m < this->nCustomerTypes; m++) {
                problem.costs[i][j][m] = new double[nTimeSteps];
                solution[i][j][m] = new int[nTimeSteps];
            }
        }
    }

    problem.app = new int***[nCells];
    for (int i = 0; i < this->nCells; i++) {
        problem.app[i] = new int**[nCells];
        for (int j = 0; j < this->nCells; j++) {
            problem.app[i][j] = new int*[nCustomerTypes];
            for (int m = 0; m < this->nCustomerTypes; m++) {
                problem.app[i][j][m] = new int[nTimeSteps];
            }
        }
    }

    problem.vetTaskOfPeople = new int[nCustomerTypes];
    problem.vetTaskCells = new int[nCells];
    problem.peopleTakenInCells = new int**[nCells];

    for (int i = 0; i < this->nCells; i++) {
        problem.peopleTakenInCells[i] = new int*[nCustomerTypes];
        for (int m = 0; m < this->nCustomerTypes; m++) {
            problem.peopleTakenInCells[i][m] = new int[nTimeSteps];
        }
    }





    problem.vetDivision = new double **[nCells];

    for ( int j = 0; j < nCells; j++ ){
        problem.vetDivision[j] = new double*[nCustomerTypes];
        for ( int i = 0; i < nCustomerTypes; i++ ){
            problem.vetDivision[j][i] = new double [nTimeSteps];

        }

    }

    /* problem.vetCellConteggio = new CellConteggio[nCells];


     for ( int cell = 0; cell < nCells; cell++ ){
         problem.vetCellConteggio[cell].vetDivision = new double *[nCustomerTypes];
         for ( int i = 0; i < nCustomerTypes; i++ ){
             problem.vetCellConteggio[cell].vetDivision[i] = new double [nTimeSteps];

         }




     }*/



    getline(iffN, line);
    getline(iffN, line);
    std::replace(line.begin(), line.end(), ';', ' ');
    istringstream issN(line);
    for (int m = 0; m < nCustomerTypes; m++) {
        issN >> word;
        problem.vetTaskOfPeople[m] = atoi(word.c_str());
    }

    getline(iffN, line);
    for (int m = 0; m < nCustomerTypes; m++) {
        for (int t = 0; t < nTimeSteps; t++) {
            getline(iffN, line);// linea con m e t
            for (int i = 0; i < nCells; i++) {
                getline(iffN, line);// linea della matrice c_{ij} per t ed m fissati
                istringstream issC(line);
                for (int j = 0; j < nCells; j++) {
                    issC >> word;
                    problem.costs[i][j][m][t] = atoi(word.c_str());


                }
            }
        }
    }

    getline(iffN, line);
    getline(iffN, line);
    std::replace(line.begin(), line.end(), ';', ' ');
    istringstream issA(line);
    for (int i = 0; i < nCells; i++) {
        issA >> word;
        problem.vetTaskCells[i] = atoi(word.c_str());
    }

    getline(iffN, line);
    for (int m = 0; m < nCustomerTypes; m++) {
        for (int t = 0; t < nTimeSteps; t++) {
            getline(iffN, line);
            getline(iffN, line);
            std::replace(line.begin(), line.end(), ';', ' ');
            istringstream issU(line);
            for (int i = 0; i < nCells; i++) {
                issU >> word;
                problem.peopleTakenInCells[i][m][t] = atoi(word.c_str());
            }
        }
    }
}


void Heuristic::solveFast(vector<double>& stat, int timeLimit) {

    //********************************************************* v7 ****************************************************
    /**                      
     *
     *  Le varie versioni sono ora chiamate come funzioni, ritornano il risultato migliore dalle loro iterazioni
     *  il tutto al fine di gestirle su thread distinti e convogliare la soluzione migliore al main
     *   versione "democratica" -> Consiste nel gestire un solo task per cella alla volta ciclando finchè non si ha raggiunta
     *   la soglia per tutte le celle
     *   vantaggi: nella v6 vengono risolti prima (e quindi in modo solitamente migliore) i task delle celle precedenti nel vettore,
     *   mentre con questa versione è possibile risolvere in modo equo ogni cella
     * Altri cambiamenti minori:
     * inserimento di una struttura dati per gestire i dati su una cella destinazione
     *
     * Altri cambiamenti minori:
     * bug fixes on feasibility check checking cell by cell
     * bug fixes all around
     */

    int bestTotalCost = INT32_MAX; //minimo globale inizializzato al max rappresentabile in intero 32 bit
    int *bestTypeTakenToStats = new int[nCustomerTypes]; //Quante persone ho preso per tipo globali
    int dimDestArray;
    int tmp,bestThread;
    int nThread=3; //flag dei thread
    int magicNumber;
    thread vThread[3]; //3 thread + il padre = 4 thread massimo
    problem.bestTypeTakenByThreads = new int*[4];
    problem.objFunOfThreads = new int[4];
    int predefinedTaskCount=problem.vetTaskOfPeople[0]; //Soglia dei task

    //**********************************Inizializzazione varie variabili***********************************************

    clock_t tStart = clock();

    for (int i = 0; i < nCells; i++) {
        for (int j = 0; j < nCells; j++) {
            for (int m = 0; m < nCustomerTypes; m++) {
                for (int tmp = 0; tmp < nTimeSteps; tmp++) {
                    solution[i][j][m][tmp] = 0;
                }
            }
        }
    }
    //Calcolo soglia
    for( tmp=1;tmp<nCustomerTypes;tmp++){
        predefinedTaskCount+=problem.vetTaskOfPeople[tmp];
    }
    //predefinedTaskCount*=3.5; //***Totalmente da esplicare @Francesco*** 3.5 - A.K.A. ***************MAGIC NUMBER****************
    //Calcolo dimensione destArray, copia di vetTaskCells e della cella correlata, allocazione copia di peopletakenInCell
    for(tmp=dimDestArray=0;tmp<nCells;tmp++){
        if(problem.vetTaskCells[tmp]>0){
            dimDestArray++;
        }
    }
    //Threads --> 3 al democratico, il padre al cellbycell
    for(tmp=0,magicNumber=3.5;tmp<nThread;tmp++){
        problem.bestTypeTakenByThreads[tmp]=new int[nCustomerTypes];
        problem.objFunOfThreads[tmp]=INT32_MAX;
        //vThread[tmp]=std::thread(&Heuristic::democratic,this,tStart,dimDestArray,tmp,predefinedTaskCount*magicNumber,0);
        //magicNumber++;
        if(tmp%2==0){ //"uno a te uno a me"
            vThread[tmp]=std::thread(&Heuristic::democratic,this,tStart,dimDestArray,tmp,predefinedTaskCount,tmp);
        }
        else{
            vThread[tmp]=std::thread(&Heuristic::cellByCell,this,tStart,dimDestArray,tmp,predefinedTaskCount,1);
        }
    }
    //Anche il padre è esso stesso un thread, quindi al lavoro!
    problem.bestTypeTakenByThreads[tmp] = new int[nCustomerTypes]; //tmp=4 padre
    problem.objFunOfThreads[tmp]=INT32_MAX;
    //democratic(tStart,dimDestArray,tmp,predefinedTaskCount*magicNumber,0);
    cellByCell(tStart,dimDestArray,tmp,predefinedTaskCount,0);
    for(tmp=0;tmp<nThread;tmp++){
        vThread[tmp].join();
    }
    for(tmp=bestThread=0;tmp<=nThread;tmp++){
        if(bestTotalCost>problem.objFunOfThreads[tmp]){
            bestThread=tmp;
            bestTotalCost=problem.objFunOfThreads[tmp];
        }
    }

    stat.push_back(bestTotalCost);
    stat.push_back((double)(clock() - tStart) / CLOCKS_PER_SEC);
    for(tmp=0;tmp<nCustomerTypes;tmp++) { //Pusho quante persone ho preso (Me lo chiede il testo)
        stat.push_back(problem.bestTypeTakenByThreads[bestThread][tmp]);
    }
    if(dimDestArray>0 && bestTotalCost == INT32_MAX){
        printf("EPIC FAIL - seppuku advertised");
    }
    hasSolution=true;
}

void Heuristic::cellByCell(double tStart,int dimDestArray,int tId,int predefinedTaskCount,int cicloN){
    //Declaration & initialization
    int tmp,tmp2,tmp3; //temporanee

    double clockCheck = 0; //controllo dei sencondi rimasti
    double *vetMinDivisione = new double[ nCustomerTypes * nCells * nTimeSteps ];
    int *vetUser = new int[ nCustomerTypes * nCells * nTimeSteps ];
    int *vetPeriod = new int[ nCustomerTypes * nCells * nTimeSteps ];
    int *vetCell = new int[ nCustomerTypes * nCells * nTimeSteps ];
    int generalIndex=0;

    int *typeTaken = new int[nCustomerTypes]; //Definisce le persone papabili ad ogni soluzione
    int *bestSolByType = new int[nCustomerTypes]; //Best solution per completare i task rimanenti
    int *typeTakenToStats = new int[nCustomerTypes]; //Quante persone ho preso per tipo (lo chiede l'output)
    int bestCost,objFun;

    CellDestData* cellUtility;
    int ***pplTakenInCellCopy = new int**[nCells]; //Ogni generazione di una soluzione avrà bisogno del suo pool di persone per essere feasible
    int *vetTaskCellsCopy = new int[dimDestArray];
    int *cellDestArray = new int[dimDestArray]; //Vettore celle destinazione che verrà utilizzato per creare pool diversi di partenza


    for(tmp=tmp2=0;tmp<nCells;tmp++){
        if(problem.vetTaskCells[tmp]>0){ //Prooning_4
            cellDestArray[tmp2]=tmp;
            vetTaskCellsCopy[tmp2++]=problem.vetTaskCells[tmp];
        }
        pplTakenInCellCopy[tmp]=new int*[nCustomerTypes];
        for(tmp3=0;tmp3<nCustomerTypes;tmp3++){
            pplTakenInCellCopy[tmp][tmp3]=new int[nTimeSteps];
        }
    }
    //Allocazione memoria per struct
    tmp=1;
    cellUtility = new CellDestData[tmp]; //Alloco il vettore di struct della dimensione giusta
    cellUtility[tmp].vetMinDivisione = new double[ nCustomerTypes * nCells * nTimeSteps ];
    cellUtility[tmp].vetUser = new int[ nCustomerTypes * nCells * nTimeSteps ];
    cellUtility[tmp].vetPeriod = new int[ nCustomerTypes * nCells * nTimeSteps ];
    cellUtility[tmp].vetCell = new int[ nCustomerTypes * nCells * nTimeSteps ];

    //******************************************************ALGORITMO*************************************************
    while(clockCheck < 4.95){ // Soglia tempo
        objFun=0;
        if(cicloN++){//****** If not the first cycle --> randomize
            //std::random_shuffle implemented for our problem
            for(int rand_dest,idx_rand = dimDestArray-1 ; idx_rand > 0 ; idx_rand--){
                rand_dest=std::rand()%(idx_rand);
                if(rand_dest!=idx_rand){
                    //swamp cell
                    tmp=cellDestArray[rand_dest];
                    cellDestArray[rand_dest]=cellDestArray[idx_rand];
                    cellDestArray[idx_rand]=tmp;
                    //swap task of cell
                    tmp=vetTaskCellsCopy[rand_dest];
                    vetTaskCellsCopy[rand_dest]=vetTaskCellsCopy[idx_rand];
                    vetTaskCellsCopy[idx_rand]=tmp;
                }
            }
        }
        for(tmp=0;tmp<nCells;tmp++){ //Reset del vettore delle persone prese a quello di partenza
            for(tmp2=0;tmp2<nCustomerTypes;tmp2++){ //necessario per mantenere le soluzioni feasible tra loro
                typeTakenToStats[tmp2]=0; //Reset statistiche della soluzione attuale (stats persone per tipo prese)
                for(tmp3=0;tmp3<nTimeSteps;tmp3++){
                    pplTakenInCellCopy[tmp][tmp2][tmp3]=problem.peopleTakenInCells[tmp][tmp2][tmp3];
                }
            }
        }
        for ( int index=0; index < dimDestArray; index++ ) {
            int cellDest = cellDestArray[index];
            int numberOfTask = vetTaskCellsCopy[index];
            if (objFun < problem.objFunOfThreads[tId]) { //Prooning_3
                for (int j = 0; j < nCustomerTypes; j++) {
                    typeTaken[j] = 0; //persone prese per tipo nella cella valutata
                }
                for (int period = 0; period < nTimeSteps; period++) {
                    for (int cellSource = 0; cellSource < nCells; cellSource++) {
                        if (cellDest != cellSource) {
                            for (int userA = 0; userA < nCustomerTypes; userA++) {
                                if (pplTakenInCellCopy[cellSource][userA][period] > 0) { //Pruning_2
                                    vetMinDivisione[generalIndex] = problem.costs[cellSource][cellDest][userA][period]/problem.vetTaskOfPeople[userA]; //Costo/n°di task che può fare
                                    vetUser[generalIndex] = userA; //qua prima o poi ci si piazza na struct
                                    vetPeriod[generalIndex] = period;
                                    vetCell[generalIndex] = cellSource;
                                    generalIndex++;
                                    typeTaken[userA] += pplTakenInCellCopy[cellSource][userA][period]; //persone papabili per tipo data una cella di destinazione
                                }
                            }
                        }
                    }
                }
                quickSort(vetMinDivisione, vetUser, vetPeriod, vetCell, 0, generalIndex - 1);
                int currentTask = 0; //task attualmente presi
                int i = 0; //indice dei vari vettori
                int numberOfPeopleTaken, taskToBeDone;
                //***** Inserisco in modo greedy fino ad una certa soglia *****
                while (currentTask < numberOfTask - predefinedTaskCount) {
                    taskToBeDone = (numberOfTask - currentTask);
                    numberOfPeopleTaken = 0;
                    while ((numberOfPeopleTaken + 1) * problem.vetTaskOfPeople[vetUser[i]] < taskToBeDone &&
                           pplTakenInCellCopy[vetCell[i]][vetUser[i]][vetPeriod[i]] > 0) {
                        pplTakenInCellCopy[vetCell[i]][vetUser[i]][vetPeriod[i]]--; //Questo mi permette di generarle già feasible
                        numberOfPeopleTaken++;
                        if (((numberOfPeopleTaken) * problem.vetTaskOfPeople[vetUser[i]] + currentTask) > (numberOfTask - predefinedTaskCount)) {
                            break; //Se supero la soglia mi fermo (interrompo il "cojo cojo"
                        }
                    }
                    currentTask += numberOfPeopleTaken * problem.vetTaskOfPeople[vetUser[i]];
                    objFun += numberOfPeopleTaken * problem.costs[vetCell[i]][cellDest][vetUser[i]][vetPeriod[i]];
                    typeTaken[vetUser[i]] -= numberOfPeopleTaken;
                    typeTakenToStats[vetUser[i]] += numberOfPeopleTaken; //Stats per l'output
                    i++;
                }
                //In realtà con quel greedy qui sopra noi abbiamo preso alla "cojo cojo" tutto il prendibile prima della soglia
                // quindi il baseIndex del powerset potrebbe partire dalla i attuale ma con un minimo controllo possiamo raccattare
                // anche le persone non ancora prese se la soglia fosse stata superata durante il while annidato
                if (i > 0) { //Optimizer_1
                    if (pplTakenInCellCopy[vetCell[i - 1]][vetUser[i - 1]][vetPeriod[i - 1]] > 0) {
                        i--;
                    }
                }
                bestCost = INT32_MAX; // massimo numero rappresentabile in int (32_bit) così siam sicuri che i task sicuramente li finisce
                for (int j = 0; j <nCustomerTypes; j++) { //massimo numero di persone del tipo j che servono per completare i task rimanenti
                    tmp = ((numberOfTask - currentTask) / problem.vetTaskOfPeople[j]) + 1;
                    if (typeTaken[j] > tmp) { //Non prendo tutte le persone papabili ma solo il massimo di quelle che realmente mi potrebbero servire
                        typeTaken[j] = tmp; //Per il caso più sfigato ovvero tutte di un tipo
                    }
                }
                //Combinazioni per elaborare i task rimanenti (post greedy)
                powerset(&numberOfTask, currentTask, vetUser, vetPeriod, vetCell, i, bestSolByType, typeTaken, &bestCost, 0, cellDest,pplTakenInCellCopy);
                if(bestCost==INT32_MAX){
                    objFun = INT32_MAX; //nota: il bestTotalCost globale è impostato inizialmente a quel valore
                    index = dimDestArray; //Prooning_5 feasibility
                }
                else {
                    // Ora dovrei sapere quante persone e di quale tipo dovrei prendere in bestSolByType (per l'esattezza in typeTaken[j] - bestSolByType[j])
                    for (int k = i, j = 0; j < nCustomerTypes; j++, k = i) { //quindi aggiorno con la soluzione
                        bestSolByType[j] = typeTaken[j] - bestSolByType[j];
                        while (bestSolByType[j] != 0) {
                            while (vetUser[k] != j ||
                                   pplTakenInCellCopy[vetCell[k]][vetUser[k]][vetPeriod[k]] <= 0) {
                                k++; //k tiene il baseIndex e viene resettato ad ogni cambio di tipo come nel powerset
                            }
                            pplTakenInCellCopy[vetCell[k]][vetUser[k]][vetPeriod[k]]--;
                            objFun += problem.costs[vetCell[k]][cellDest][vetUser[k]][vetPeriod[k]];
                            currentTask += problem.vetTaskOfPeople[vetUser[k]];
                            bestSolByType[j]--;
                            typeTakenToStats[j]++; //stats per l'output
                        }
                    }
                    if (currentTask < numberOfTask) { //Non feasible - ho finito le persone mentre cercavo di completare i task della cella
                        objFun = INT32_MAX; //nota: il bestTotalCost globale è impostato inizialmente a quel valore
                        index = dimDestArray; //Prooning_5 feasibility
                    }
                }
            }
            else{
                objFun= INT32_MAX; //Prooning_3
                index=dimDestArray;
            }
            generalIndex = 0;
        }
        if (problem.objFunOfThreads[tId] > objFun) { //Aggiorno la soluzione migliore
            problem.objFunOfThreads[tId] = objFun;
            for (tmp = 0; tmp < nCustomerTypes; tmp++) {
                problem.bestTypeTakenByThreads[tId][tmp] = typeTakenToStats[tmp];
            }
        }
        clockCheck = (double)(clock() - tStart) / CLOCKS_PER_SEC;
    }
    return;
}

void Heuristic::democratic(double tStart,int dimDestArray,int tId,int predefinedTaskCount,int cicloN){
    //Declaration & initialization
    int tmp,tmp2,tmp3; //temporanee

    double clockCheck = 0; //controllo dei sencondi rimasti

    int *typeTaken = new int[nCustomerTypes]; //Definisce le persone papabili ad ogni soluzione
    int *bestSolByType = new int[nCustomerTypes]; //Best solution per completare i task rimanenti
    int *typeTakenToStats = new int[nCustomerTypes]; //Quante persone ho preso per tipo (lo chiede l'output)
    int bestCost,objFun;

    CellDestData* cellUtility;
    int ***pplTakenInCellCopy = new int**[nCells]; //Ogni generazione di una soluzione avrà bisogno del suo pool di persone per essere feasible
    int *typeTakenRemaining = new int[nCustomerTypes];
    int *vetTaskCellsCopy = new int[dimDestArray];
    int *cellDestArray = new int[dimDestArray]; //Vettore celle destinazione che verrà utilizzato per creare pool diversi di partenza

    cellUtility = new CellDestData[dimDestArray]; //Alloco il vettore di struct della dimensione giusta

    for(tmp=tmp2=0;tmp<nCells;tmp++){
        if(problem.vetTaskCells[tmp]>0){ //Prooning_4
            cellDestArray[tmp2]=tmp;
            vetTaskCellsCopy[tmp2++]=problem.vetTaskCells[tmp];
        }
        pplTakenInCellCopy[tmp]=new int*[nCustomerTypes];
        for(tmp3=0;tmp3<nCustomerTypes;tmp3++){
            pplTakenInCellCopy[tmp][tmp3]=new int[nTimeSteps];
        }
    }
    for(tmp=0;tmp<dimDestArray;tmp++){
        //Allocazione memoria per struct
        cellUtility[tmp].vetMinDivisione = new double[ nCustomerTypes * nCells * nTimeSteps ];
        cellUtility[tmp].vetUser = new int[ nCustomerTypes * nCells * nTimeSteps ];
        cellUtility[tmp].vetPeriod = new int[ nCustomerTypes * nCells * nTimeSteps ];
        cellUtility[tmp].vetCell = new int[ nCustomerTypes * nCells * nTimeSteps ];
    }
    //******************************************************ALGORITMO*************************************************
    while(clockCheck < 4.95){ // Soglia tempo
        objFun=0;
        if(cicloN++){//****** If not daddy --> randomize
            //std::random_shuffle implemented for our problem
            for(int rand_dest,idx_rand = dimDestArray-1 ; idx_rand > 0 ; idx_rand--){
                rand_dest=std::rand()%(idx_rand);
                if(rand_dest!=idx_rand){
                    //swamp cell
                    tmp=cellDestArray[rand_dest];
                    cellDestArray[rand_dest]=cellDestArray[idx_rand];
                    cellDestArray[idx_rand]=tmp;
                    //swap task of cell
                    tmp=vetTaskCellsCopy[rand_dest];
                    vetTaskCellsCopy[rand_dest]=vetTaskCellsCopy[idx_rand];
                    vetTaskCellsCopy[idx_rand]=tmp;
                }
            }
        }
        for(tmp=0;tmp<nCells;tmp++){ //Reset del vettore delle persone prese a quello di partenza
            for(tmp2=0;tmp2<nCustomerTypes;tmp2++){ //necessario per mantenere le soluzioni feasible tra loro
                typeTakenToStats[tmp2]=0; //Reset statistiche della soluzione attuale (stats persone per tipo prese)
                typeTaken[tmp2]=0;//reset persone della soluzione
                for(tmp3=0;tmp3<nTimeSteps;tmp3++){
                    pplTakenInCellCopy[tmp][tmp2][tmp3]=problem.peopleTakenInCells[tmp][tmp2][tmp3];
                }
            }
        }
        //************************************************SOLUZIONE DEMOCRATICA****************************************
        //                              INIZIALIZZAZIONE DELLE CELLE (STRUCT) PER LA SOLUZIONE CORRENTE
        //Generazione vettore di nCells strutture contente i dati utili per ogni cella & reset statistiche della soluzione corrente
        for ( int index=0; index < dimDestArray; index++ ) {
            cellUtility[index].cellDest = cellDestArray[index];
            cellUtility[index].numberOfTask = vetTaskCellsCopy[index];
            cellUtility[index].currentTask = 0;
            cellUtility[index].startingPoint = 0;
            cellUtility[index].flag = 0;
            cellUtility[index].size = 0; //quante persone ho nel vettore ordinato
            for (int period = 0; period < nTimeSteps; period++) {
                for (int cellSource = 0; cellSource < nCells; cellSource++) {
                    if (cellUtility[index].cellDest != cellSource) {
                        for (int userA = 0; userA < nCustomerTypes; userA++) {
                            if (pplTakenInCellCopy[cellSource][userA][period] > 0) { //Pruning_2
                                cellUtility[index].vetMinDivisione[cellUtility[index].size] = problem.costs[cellSource][cellUtility[index].cellDest][userA][period] /
                                                                                              problem.vetTaskOfPeople[userA]; //Costo/n°di task che può fare
                                cellUtility[index].vetUser[cellUtility[index].size] = userA;
                                cellUtility[index].vetPeriod[cellUtility[index].size] = period;
                                cellUtility[index].vetCell[cellUtility[index].size] = cellSource;
                                cellUtility[index].size++;
                                typeTaken[userA] += pplTakenInCellCopy[cellSource][userA][period]; //persone papabili per tipo data una cella di destinazione
                            }
                        }
                    }
                }
            }
            quickSort(cellUtility[index].vetMinDivisione, cellUtility[index].vetUser, cellUtility[index].vetPeriod, cellUtility[index].vetCell, 0, cellUtility[index].size - 1);
        }
        //                                                              **PARTE GREEDY**
        //Ciclo sulle celle risolvendo un task alla volta
        //Condizione di terminazione: soglie raggiunte = nCells
        int countSolved = 0;
        int dem_index;
        int i = 0;
        while(countSolved < dimDestArray){
            for(dem_index = 0; dem_index < dimDestArray ; dem_index++){//Ciclo sulle celle destinazione
                if(cellUtility[dem_index].flag == 0){ //skippo se ho già risolto
                    if(cellUtility[dem_index].currentTask < cellUtility[dem_index].numberOfTask - predefinedTaskCount){// Risolvo un task alla volta
                        i=cellUtility[dem_index].startingPoint;
                        while(pplTakenInCellCopy[cellUtility[dem_index].vetCell[i]][cellUtility[dem_index].vetUser[i]][cellUtility[dem_index].vetPeriod[i]] < 0){ //Trovo la prima persona disponibiòe
                            i++;
                        }
                        pplTakenInCellCopy[cellUtility[dem_index].vetCell[i]][cellUtility[dem_index].vetUser[i]][cellUtility[dem_index].vetPeriod[i]]--; //Tolgo la persona dall'elenco
                        cellUtility[dem_index].currentTask += problem.vetTaskOfPeople[cellUtility[dem_index].vetUser[i]]; //Aggiungo i task che fa quell'utente
                        cellUtility[dem_index].startingPoint = i; //Aggiorno il punto di partenza del powerset
                        objFun += problem.costs[cellUtility[dem_index].vetCell[i]][cellUtility[dem_index].cellDest][cellUtility[dem_index].vetUser[i]][cellUtility[dem_index].vetPeriod[i]];
                        typeTaken[cellUtility[dem_index].vetUser[i]]--; //Tolgo dall'elenco delle persone papabili per tipo
                        typeTakenToStats[cellUtility[dem_index].vetUser[i]]++;
                    }
                    else{
                        countSolved++;
                        cellUtility[dem_index].flag = 1;
                    }
                }
            }
            dem_index = 0;
        }
        //                                                **Task sotto soglia --> powerset**
        for ( int index=0; index < dimDestArray; index++ ) { //Scorro sulle celle già risolte fino alla soglia
            if (objFun < problem.objFunOfThreads[tId]) { //Prooning_3
                if(cellUtility[index].numberOfTask-cellUtility[index].currentTask > 0) { //Prooning
                    bestCost = INT32_MAX; // massimo numero rappresentabile in int (32_bit) così siam sicuri che i task sicuramente li finisce
                    for (tmp = 0; tmp <
                                  nCustomerTypes; tmp++) { //massimo numero di persone del tipo j che servono per completare i task rimanenti
                        tmp2 = ((cellUtility[index].numberOfTask - cellUtility[index].currentTask) /
                                problem.vetTaskOfPeople[tmp]) + 1;
                        if (typeTaken[tmp] >
                            tmp2) { //Non prendo tutte le persone papabili ma solo il massimo di quelle che realmente mi potrebbero servire
                            typeTakenRemaining[tmp] = tmp2; //Per il caso più sfigato ovvero tutte di un tipo
                        } else {
                            typeTakenRemaining[tmp] = typeTaken[tmp]; //Vado a modificare i minimi da prendere per questa cella
                        }
                    }
                    //Combinazioni per elaborare i task rimanenti (post greedy)
                    powerset(&cellUtility[index].numberOfTask, cellUtility[index].currentTask,
                             cellUtility[index].vetUser, cellUtility[index].vetPeriod,
                             cellUtility[index].vetCell, cellUtility[index].startingPoint, bestSolByType,
                             typeTakenRemaining, &bestCost, 0,
                             cellUtility[index].cellDest, pplTakenInCellCopy);
                    // Ora dovrei sapere quante persone e di quale tipo dovrei prendere in bestSolByType (per l'esattezza in typeTaken[j] - bestSolByType[j])
                    if (bestCost == INT32_MAX) {
                        objFun = bestCost; //Se non ho trovato una soluzione feasible nel powerset questa soluzione non va bene
                        index = dimDestArray;
                    } else {
                        for (int k = cellUtility[index].startingPoint, j = 0; j <
                                                                              nCustomerTypes; j++, k = cellUtility[index].startingPoint) { //quindi aggiorno con la soluzione
                            bestSolByType[j] = typeTakenRemaining[j] - bestSolByType[j];
                            while (bestSolByType[j] != 0) {
                                while (cellUtility[index].vetUser[k] != j ||
                                       pplTakenInCellCopy[cellUtility[index].vetCell[k]][cellUtility[index].vetUser[k]][cellUtility[index].vetPeriod[k]] <=
                                       0) {
                                    k++; //k tiene il baseIndex e viene resettato ad ogni cambio di tipo come nel powerset
                                }
                                pplTakenInCellCopy[cellUtility[index].vetCell[k]][cellUtility[index].vetUser[k]][cellUtility[index].vetPeriod[k]]--;
                                objFun += problem.costs[cellUtility[index].vetCell[k]][cellUtility[index].cellDest][cellUtility[index].vetUser[k]][cellUtility[index].vetPeriod[k]];
                                cellUtility[index].currentTask += problem.vetTaskOfPeople[cellUtility[index].vetUser[k]];
                                bestSolByType[j]--;
                                typeTaken[j]--;
                                typeTakenToStats[j]++; //stats per l'output
                            }
                        }
                        if (cellUtility[index].currentTask <
                            cellUtility[index].numberOfTask) { //Non feasible - ho finito le persone mentre cercavo di completare i task della cella
                            objFun = INT32_MAX; //nota: il bestTotalCost globale è impostato inizialmente a quel valore
                            index = dimDestArray; //Prooning_5 feasibility
                        }
                    }
                }
            }
            else{
                objFun= INT32_MAX; //Prooning_3
                index=dimDestArray;
            }
        }
        if (problem.objFunOfThreads[tId] > objFun) { //Aggiorno la soluzione migliore
            problem.objFunOfThreads[tId] = objFun;
            for (tmp = 0; tmp < nCustomerTypes; tmp++) {
                problem.bestTypeTakenByThreads[tId][tmp] = typeTakenToStats[tmp];
            }
        }
        clockCheck = (double)(clock() - tStart) / CLOCKS_PER_SEC;
    }
    return;
}

void Heuristic::powerset(int *numberOfTask,int currentTask,int* vetUser,int* vetPeriod,int* vetCell,int baseIndex,int *bestSolByType,
                         int* typeTaken,int *bestCost,int currCost,int cellDest,int ***pplTakenInCellCopy){
    int k,j,i;
    if (currCost > *bestCost ){ //Prooning_0
        return;
    }
    else if(currentTask >= *numberOfTask){ //terminazione
        if(*bestCost==currCost){ //Se ho preso meno persone a pari costo ok, altrimenti return
            for(j=i=k=0;j<nCustomerTypes;j++){
                i+=bestSolByType[j];
                k+=typeTaken[j];
            }
            if(k>i){
                return;
            }
        }
        *bestCost = currCost; //Se currCost fosse stato più grande di bestCost sarebbe returnato prima
        for(j=0;j<nCustomerTypes; j++) { //aggiorno soluzione migliore
            bestSolByType[j] = typeTaken[j];
        }
        return;
    }
    for(j=0,i=baseIndex; j<nCustomerTypes; j++, i=baseIndex) { //Ad ogni cambio del tipo analizzato torno al baseIndex
        if (typeTaken[j] > 0) {
            typeTaken[j]--; //prendo una persona in più dal pool (se c'è) ed effettuo la ricorsione
            while (vetUser[i] != j || pplTakenInCellCopy[vetCell[i]][vetUser[i]][vetPeriod[i]] <= 0) {
                i++;
            }
            pplTakenInCellCopy[vetCell[i]][vetUser[i]][vetPeriod[i]]--; //Ricorsione
            powerset(numberOfTask,
                     (currentTask + problem.vetTaskOfPeople[vetUser[i]]), //New value for currTask
                     vetUser, vetPeriod, vetCell, baseIndex, bestSolByType,
                     typeTaken, bestCost,
                     (currCost + problem.costs[vetCell[i]][cellDest][vetUser[i]][vetPeriod[i]]), //New value for currCost
                     cellDest,pplTakenInCellCopy);
            pplTakenInCellCopy[vetCell[i]][vetUser[i]][vetPeriod[i]]++; //Reset post ricorsione
            typeTaken[j]++;
        }
    }
    return;
}


void Heuristic::quickSort2(int *arr, int *arr2, int left,int right) // --- Fail_1 ---
{
    int i = left, j = right;
    int tmp;
    int pivot = arr[(left + right) / 2];

    /* partition */
    while (i <= j) {
        while (arr[i] < pivot)
            i++;
        while (arr[j] > pivot)
            j--;
        if (i <= j) {
            tmp = arr[i];
            arr[i] = arr[j];
            arr[j] = tmp;

            tmp = arr2[i];
            arr2[i] = arr2[j];
            arr2[j] = tmp;

            i++;
            j--;
        }
    };

    /* recursion */
    if (left < j)
        quickSort2(arr, arr2 , left, j);
    if (i < right)
        quickSort2(arr, arr2 , i, right);
}



void Heuristic::quickSort(double *arr, int *arr2 , int *arr3 , int *arr4 , int left, int right) {
    int i = left, j = right;
    double tmp;
    int tmp2 ;
    double pivot = arr[(left + right) / 2];

    /* partition */
    while (i <= j) {
        while (arr[i] < pivot)
            i++;
        while (arr[j] > pivot)
            j--;
        if (i <= j) {
            tmp = arr[i];
            arr[i] = arr[j];
            arr[j] = tmp;

            tmp2 = arr2[i];
            arr2[i] = arr2[j];
            arr2[j] = tmp2;

            tmp2 = arr3[i];
            arr3[i] = arr3[j];
            arr3[j] = tmp2;

            tmp2 = arr4[i];
            arr4[i] = arr4[j];
            arr4[j] = tmp2;

            i++;
            j--;
        }
    };

    /* recursion */
    if (left < j)
        quickSort(arr, arr2 , arr3 , arr4 , left, j);
    if (i < right)
        quickSort(arr, arr2 , arr3 , arr4 , i, right);
}



void Heuristic::writeKPI(string path, string instanceName, vector<double> stat){
    if (!hasSolution)
        return;

    ofstream fileO(path, ios::app);
    if(!fileO.is_open())
        return;

    fileO << instanceName << ";" << stat[0] << ";" << stat[1];
    for(int i=2; i<nCustomerTypes+2; i++)
        fileO <<  ";" << stat[i];
    fileO << endl;

    fileO.close();

}

void Heuristic::writeSolution(string path) {
    if (!hasSolution)
        return;

    ofstream fileO(path);
    if(!fileO.is_open())
        return;

    fileO << this->nCells << "; " << this->nTimeSteps << "; " << this->nCustomerTypes << endl;
    for (int m = 0; m < this->nCustomerTypes; m++)
        for (int t = 0; t < this->nTimeSteps; t++)
            for (int i = 0; i < this->nCells; i++)
                for (int j = 0; j < this->nCells; j++)
                    if (solution[i][j][m][t] > 0)
                        fileO << i << ";" << j << ";" << m << ";" << t << ";" << solution[i][j][m][t] << endl;

    fileO.close();
}

eFeasibleState Heuristic::isFeasible(string path) {

    string line;
    string word;
    int nCellsN;
    int nTimeStepsN;
    int nCustomerTypesN;
    int i, j, m, t;


    ifstream iffN(path.c_str());

    if (!iffN.is_open()) {
        cout << "Impossible to open" << path << endl;
        exit(1);
    }

    getline(iffN, line);
    std::replace(line.begin(), line.end(), ';', ' ');
    istringstream iss(line);
    iss >> word; // nCells
    nCellsN = atoi(word.c_str());
    iss >> word; // nTimeSteps
    nTimeStepsN = atoi(word.c_str());
    iss >> word; // nCustomerTypes
    nCustomerTypesN = atoi(word.c_str());

    int**** solutionN = new int***[nCells];
    for (i = 0; i < nCellsN; i++) {
        solutionN[i] = new int**[nCells];
        for (j = 0; j < nCellsN; j++) {
            solutionN[i][j] = new int*[nCustomerTypes];
            for (m = 0; m < nCustomerTypesN; m++) {
                solutionN[i][j][m] = new int[nTimeSteps];
                for ( t = 0; t < nTimeStepsN; t++) {
                    solutionN[i][j][m][t] = 0;
                }
            }
        }
    }

    while (getline(iffN, line)) {
        std::replace(line.begin(), line.end(), ';', ' ');
        istringstream iss(line);
        iss >> word; // i
        i = atoi(word.c_str());
        iss >> word; // j
        j = atoi(word.c_str());
        iss >> word; // m
        m = atoi(word.c_str());
        iss >> word; // t
        t = atoi(word.c_str());
        iss >> word; // value
        solutionN[i][j][m][t] = atoi(word.c_str());
    }

    // Demand
    bool feasible = true;
    int expr = 0;
    for (int i = 0; i < nCells; i++) {
        for (int j = 0; j < nCells; j++)
            for (int m = 0; m < nCustomerTypes; m++)
                for (int t = 0; t < nTimeSteps; t++)
                    expr += problem.vetTaskOfPeople[m] * solutionN[j][i][m][t];
        if (expr < problem.vetTaskCells[i])
            feasible = false;
    }

    if (!feasible)
        return NOT_FEASIBLE_DEMAND;

    // Max Number of users
    for (i = 0; i < nCells; i++)
        for (m = 0; m < nCustomerTypes; m++)
            for (t = 0; t < nTimeSteps; t++) {
                expr = 0;
                for (j = 0; j < nCells; j++)
                    expr += solutionN[i][j][m][t];
                if (expr > problem.peopleTakenInCells[i][m][t])
                    feasible = false;
            }

    if(!feasible)
        return NOT_FEASIBLE_USERS;

    return FEASIBLE;
}

void Heuristic::getStatSolution(vector<double>& stat) {
    if (!hasSolution)
        return;

    int* tipi = new int[nCustomerTypes];
    for (int m = 0; m < nCustomerTypes; m++)
        tipi[m] = 0;

    for (int i = 0; i < nCells; i++)
        for (int j = 0; j < nCells; j++)
            for (int t = 0; t < nTimeSteps; t++)
                for (int m = 0; m < nCustomerTypes; m++)
                    if (solution[i][j][m][t] > 0)
                        tipi[m] += solution[i][j][m][t];
    for (int m = 0; m < nCustomerTypes; m++)
        stat.push_back(tipi[m]);

}
