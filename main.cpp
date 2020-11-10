#include <iostream>
#include <random>
#include <vector>
#include <algorithm>
#include <string>
#include <iomanip>

using std::cout;
using std::vector;

const double X_LOWER_EDGE = -2.;
const double X_UPPER_EDGE = 2.;
const double Y_LOWER_EDGE = -2.;
const double Y_UPPER_EDGE = 2.;
const double EPSILON_FOR_MUTATION = 0.4;
const double PROBABILITY_OF_MUTATION = 0.3;

double FitnessFunction(const double &x, const double &y) {
    return sin(x) * exp(-x * x - y * y) + 0.5;
}
// у исходной функции глоб. максимум равен 0.396653 в точке (0.653271, 0)
// так как на обл. опред. исходная функция принимает отрицательные значения,
// то в качестве FitnessFunction возьмем функцию, "поднятую" на 0,5

class Specimen {
public:
    double x;
    double y;
    double valueOfFitnessFunction;

    Specimen() = default;

    Specimen(const double &chromosome1, const double &chromosome2) : x(
            chromosome1), y(chromosome2) {
        valueOfFitnessFunction = FitnessFunction(x, y);
    }

    ~Specimen() = default;

    void UpdateValueOfFitnessFunction() {
        valueOfFitnessFunction = FitnessFunction(x, y);
    }
};

void PrintTopLine(std::stringstream &out) {
    out << std::setw(3) << " N  "
        << std::setw(9) << "x  "
        << std::setw(11) << "y  "
        << std::setw(11) << "FIT  "
        << std::setw(13) << "Max FIT  "
        << std::setw(12) << "Average FIT\n";
}

void PrintPopulation(const int &number, const vector<Specimen> &population,
                     std::stringstream &out) {
    double average = 0;
    for (const auto &specimen : population) {
        average += specimen.valueOfFitnessFunction;
    }
    average /= population.size();
    for (auto specimen = population.begin();
         specimen != population.end(); specimen++) {
        out << std::setw(2)
            << (specimen == population.begin() ? std::to_string(number)
                                               : " ");
        out << "  " << std::fixed << std::setprecision(7)
            << std::setw(10) << specimen->x << ' '
            << std::setw(10) << specimen->y << ' '
            << std::setw(10) << specimen->valueOfFitnessFunction << ' ';
        if (specimen == population.begin())
            out << std::setw(10) << specimen->valueOfFitnessFunction << ' '
                << std::setw(10) << average;
        out << '\n';
    }
}

vector<std::pair<Specimen, Specimen>> Selection(vector<Specimen> &population) {
    vector<std::pair<Specimen, Specimen>> pairsForCrossover{
            std::make_pair(population[0], population[1]),
            std::make_pair(population[0], population[2])};
    return pairsForCrossover;
}

vector<Specimen> Crossover(const vector<std::pair<Specimen, Specimen>> &pairs) {
    vector<Specimen> population(2 * pairs.size());
    size_t i = 0;
    for (auto &pair : pairs) {
        population[i] = Specimen(pair.first.x, pair.second.y);
        ++i;
        population[i] = Specimen(pair.second.x, pair.first.y);
        ++i;
    }
    return population;
}

void Mutation(vector<Specimen> &population, std::mt19937 &gen) {
    for (auto &specimen: population) {
        std::discrete_distribution<> distrib(
                {1 - PROBABILITY_OF_MUTATION, PROBABILITY_OF_MUTATION});
        if (distrib(gen)) {
            std::uniform_real_distribution<double> distForDelta(
                    -EPSILON_FOR_MUTATION / 2, EPSILON_FOR_MUTATION / 2);
            double deltaX = distForDelta(gen);
            //если вышли за обл. опред. FIT-функции
            while (specimen.x + deltaX < X_LOWER_EDGE ||
                   specimen.x > X_UPPER_EDGE) {
                deltaX = distForDelta(gen);
            }
            specimen.x += deltaX;

            double deltaY = distForDelta(gen);
            while (specimen.y + deltaY < Y_LOWER_EDGE ||
                   specimen.y > Y_UPPER_EDGE) {
                deltaY = distForDelta(gen);
            }
            specimen.y += deltaY;
            specimen.UpdateValueOfFitnessFunction();
        }
    }
}

void Reduction(vector<Specimen> &population) {
    //Удаление старых особей
    size_t sizeOfPopulation = population.size();
    population.erase(population.begin(),
                     population.begin() + sizeOfPopulation / 2);
    //сортировка для правильного вывода и дальнейшей селекции
    std::sort(std::begin(population), std::end(population),
              [](const Specimen &a, const Specimen &b) {
                  return a.valueOfFitnessFunction > b.valueOfFitnessFunction;
              });
}

Specimen ClassicGeneticAlgorithm() {
    std::stringstream table1, table2; //для 10 итерация и
                                    // для 100 итераций с шагом 10
    PrintTopLine(table1);
    PrintTopLine(table2);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> Xdist(X_LOWER_EDGE, X_UPPER_EDGE);
    std::uniform_real_distribution<double> Ydist(Y_LOWER_EDGE, Y_UPPER_EDGE);
    vector<Specimen> population(4);
    for (auto &specimen : population) {
        specimen = Specimen(Xdist(gen), Ydist(gen));
    }
    std::sort(std::begin(population), std::end(population),
              [](const Specimen &a, const Specimen &b) {
                  return a.valueOfFitnessFunction > b.valueOfFitnessFunction;
              });
    PrintPopulation(0, population, table1);
    PrintPopulation(0, population, table2);
    for (int number = 1; number < 101; ++number) {
        //Селекция - выбор хромосом, которые будут принимать участие в скрещивании
        //Вероятность выбора хромосом пропорциональна значению FIT-функции
        vector<std::pair<Specimen, Specimen>> pairsForCrossover = Selection(
                population);
        //кроссовер - скрещивание
        vector<Specimen> newPopulation = Crossover(pairsForCrossover);
        for (const auto &newSpecimen : newPopulation) {
            population.push_back(newSpecimen);
        }
        //мутация с вероятностью 25%
        Mutation(population, gen);
        //Редукция
        Reduction(population);
        if (number < 11) PrintPopulation(number, population, table1);
        if (number % 10 == 0) PrintPopulation(number, population, table2);
    }
    cout << "First 10 iteration :\n" << table1.str() << '\n';
    cout << "Every 10th iteration :\n" << table2.str() << '\n';
    return population[0];
}

int main() {
    cout << "\tinitial function : sin(x) * exp(-x * x - y * y)\n";
    cout << "\tfitness function : sin(x) * exp(-x * x - y * y) + 0.5\n";
    Specimen bestSpecimen = ClassicGeneticAlgorithm();
    cout << "\tmax of fitness function : "
         << bestSpecimen.valueOfFitnessFunction << '\n';
    cout << "\tmax of initial function : "
         << bestSpecimen.valueOfFitnessFunction - 0.5;
    cout << "   at point = ( " << bestSpecimen.x << ", " << bestSpecimen.y
         << " )\n";
    cout << "\treal max of initial function : " << 0.396653;
    cout << "   at point = ( " << 0.653271 << ", " << 0 << " )\n";
}
