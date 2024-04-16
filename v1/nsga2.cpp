#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>

using namespace std;


const int numObjectives = 2;                            // Número de objtivos:
const int populationSize = 10;                         // Tamanho da população:
const int numGenerations = 10;                         // Por exemplo, 100 gerações


struct Solution {   
    vector<double> objectives;                          // Valores dos objetivos
    double fitness;                                     // valores do fitness
};

// Função para inicializar a população com soluções aleatórias
vector<Solution> initializePopulation() {
    vector<Solution> population(populationSize);

    for (int i = 0; i < populationSize; i++) {
        population[i].objectives.resize(numObjectives);

        for (int j = 0; j < numObjectives; j++) {
            population[i].objectives[j] = (double)rand() / RAND_MAX; //Gerar valores aleatórios para os objetivos
        };

        population[i].fitness = 0.0;    // Inicializando o fitness
    };

    return population;
};

// Função para avaliar a dominância entre duas soluções
bool dominates(const Solution& a, const Solution& b) {
    bool betterInAtLeastOneObjective = false;

    for (int i = 0; i < numObjectives; i++) {
        if (a.objectives[i] < b.objectives[i]) {
            betterInAtLeastOneObjective = true;
        }
        else if (a.objectives[i] > b.objectives[i]) {
            return false;                               // Não domina se for pior em pelo menos um objetivo
        };
    };

    return betterInAtLeastOneObjective;
};

// Função para calcular o ranking de dominância para a população
vector<int> calculateDominanceRanking(const vector<Solution>& population) {
    vector<int> ranking(populationSize, 0);

    for (int i = 0; i < populationSize; i++) {
        for (int j = 0; j < populationSize; j++) {
            if (dominates(population[i], population[j])) {
                ranking[j]++;
            };
        };
    };
    
    return ranking;
};


// Função para encontrar o índice de um elemento em um vetor
template <typename T>
int findIndex(const std::vector<T>& vec, const T& value) {
    auto it = std::find(vec.begin(), vec.end(), value);
    if (it != vec.end()) {
        return std::distance(vec.begin(), it);
    };
    return -1; // Se o elemento não for encontrado, retorna -1
};


// Função para seleção por torneio
Solution tournamentSelection(const vector<Solution>& population, const vector<int>& dominanceRanking) {
    const int tournamentSize = 2; // Tamanho do torneio
    Solution bestSolution = population[rand() % populationSize]; // Escolhe um indivíduo aleatório para iniciar

    for (int i = 1; i < tournamentSize; ++i) {
        Solution challenger = population[rand() % populationSize]; // Escolhe outro indivíduo aleatório
        
        // Compara o desempenho do desafiante com o melhor indivíduo atual
        if (dominanceRanking[findIndex(population, challenger)] < dominanceRanking[findIndex(population, bestSolution)]) {
            bestSolution = challenger; // Atualiza o melhor indivíduo se o desafiante for melhor
        };
    };

    return bestSolution;
};


// Função para realizar crossover uniforme
Solution uniformCrossover(const Solution& parent1, const Solution& parent2) {
    Solution offspring;
    offspring.objectives.resize(numObjectives);

    for (int i = 0; i < numObjectives; i++) {
        if (rand() % 2 == 0) {
            offspring.objectives[i] = parent1.objectives[i];
        }
        else {
            offspring.objectives[i] = parent2.objectives[i];
        };
    };

    return offspring;
}


// Função para realizar mutação
void mutate(Solution& solution, double mutationRate) {
    for (int i = 0; i < numObjectives; i++) {
        if (rand() / static_cast<double>(RAND_MAX) < mutationRate) {
            solution.objectives[i] += ((double)rand() / RAND_MAX - .5) * .1;    // Mutação simples
        };
    };
};


// Função para calcular o fitness de cada solução
void calculateFitness(vector<Solution>& population, const vector<int>& dominanceRanking) {
    for (int i = 0; i < populationSize; i++) {
        int dominanceCount = 0;

        for (int j = 0; j < populationSize; j++) {
            if (dominanceRanking[j] < dominanceRanking[i]) {
                dominanceCount++;
            };
        };

        population[i].fitness = dominanceCount;
    };
};

template <typename T>
int findIndex(const std::vector<T>& vec, const T& value) {
    auto it = std::find(vec.begin(), vec.end(), value);
    if (it != vec.end()) {
        return std::distance(vec.begin(), it);
    }
    return -1; // Se o elemento não for encontrado, retorna -1
}

// NSGA-II
void nsga2() {
    vector<Solution> population = initializePopulation();       // Inicialização da população

    for (int generation = 0; generation < numGenerations; generation++) {       // Loop principal (número fixo de gerações)
        vector<int> dominanceRanking = calculateDominanceRanking(population);           // Avaliação de dominância e ranking


        // Seleção e crossover
        vector<Solution> offspringPopulation;
        for (int i = 0; i < populationSize; i++) {
            Solution parent1    = tournamentSelection(population, dominanceRanking);
            Solution parent2    = tournamentSelection(population, dominanceRanking);
            Solution offspring  = uniformCrossover(parent1, parent2); 

            // mutação
            mutate(offspring, .1);  // Taxa de mutação = 10%

            offspringPopulation.push_back(offspring);
        };


        // Avaliação de fitness
        calculateFitness(offspringPopulation, dominanceRanking);

        // Seleção da próxima geração
        vector<Solution> nextGeneration;
        for (int i = 0; i < populationSize; i++) {
            Solution bestSolution = tournamentSelection(population, dominanceRanking);
            nextGeneration.push_back(bestSolution);
        };


        population = nextGeneration;
    };

};


int main() {
    srand(time(nullptr)); // Inicialização da semente aleatória

    nsga2(); // Chamar a função principal do NSGA-II

    return 0;
}