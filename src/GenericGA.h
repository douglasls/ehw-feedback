#ifndef _GENERIC_GA_H
#define _GENERIC_GA_H

#include <vector>
#include <functional>
#include <math.h>
#include "random.h"
#include "Utils.h"
#include <chrono>

using namespace std;
using namespace std::chrono;

template <typename T>
struct Evaluated {
	double score;
	T value;
};

template <typename T>
struct GAState {
	int generation;
	std::vector<T> population;
};

template <typename T>
Evaluated<T> __attribute__((pure)) minimumEvaluated(std::vector<Evaluated<T>> v) {
	return fold1([](Evaluated<T> best, Evaluated<T> current) {
		if (current.score <= best.score) {
			return current;
		}
		return best;
	}, v);
}

template <typename T>
Evaluated<T> __attribute__((pure)) makeEvaluated(T value, double score) {
	Evaluated<T> res;
	res.value = value;
	res.score = score;
	return res;
}

template <typename T, typename F, typename G>
std::function<RNGFUNC(std::vector<Evaluated<T>>)(std::vector<Evaluated<T>>)>
	lambdaPlusN
			( F fitness
			, G mutation
			, int n
			) {

	static_assert(std::is_convertible<F, std::function<double(T)>> ::value,
			"lamdaPlusN's fitness function must be of type T -> double");
	static_assert(std::is_convertible<G, std::function<RNGFUNC(T)(T)>> ::value,
			"lambdaPlusN's mutation function must be of type T -> RNGFUNC(T)");

	return [=](std::vector<Evaluated<T>> population) {
		
		high_resolution_clock::time_point t1 = high_resolution_clock::now(); // INICIO TEMPO DE GERAR INDIVIDUOS
		auto result = mapM(mutation, replicate(n, population[0].value));
		high_resolution_clock::time_point t2 = high_resolution_clock::now(); //FIM
		auto TempIndiv = duration_cast<microseconds>( t2 - t1 ).count();
		std::cout << "Tempo para gerar individuo: " << TempIndiv << std::endl;
		
		return bind(result,
				[=](std::vector<T> newIndividuals) {
		    newIndividuals.insert(newIndividuals.begin(), population[0].value);
			
			high_resolution_clock::time_point t3 = high_resolution_clock::now(); // INICIO TEMPO DE AVALIAR OS INDIVIDUOS
			auto fitnesses = map(fitness, newIndividuals);
			high_resolution_clock::time_point t4 = high_resolution_clock::now(); // FIM
			auto TempAvalia = duration_cast<microseconds>( t4 - t3 ).count();
			std::cout << "Tempo para avaliar individuo: " << TempAvalia << std::endl;

			auto evaluated = zipWith(makeEvaluated<T>, newIndividuals, fitnesses);
            return pure(std::vector<Evaluated<T>>{ minimumEvaluated(evaluated) });
		});
	};
}

template <typename T, typename F>
std::function<RNGFUNC(GAState<T>)(GAState<T>)> __attribute__((pure)) makeGAFunction
    (F strategy) {

	static_assert(std::is_convertible<F, std::function<RNGFUNC(std::vector<T>)(std::vector<T>)>> ::value,
			"makeGAFunction's strategy function must be of type std::vector<T> -> RNGFUNC(std::vector<T>)");

	return [strategy](GAState<T> state) {
		return bind
				( strategy(state.population)
				, [=](std::vector<T> newPopulation) {
			GAState<T> newState;
			newState.population = newPopulation;
			newState.generation = state.generation + 1;
			return pure(newState);
		});
	};
}

#endif
