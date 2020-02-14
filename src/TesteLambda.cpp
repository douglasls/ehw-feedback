//============================================================================
// Name        : TesteLambda.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <iomanip>
#include <fstream>
#include <functional>
#include <vector>
#include <bitset>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <algorithm>
#include <tuple>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <stdint.h>
#include <stdio.h>
#include <cassert>
#include "GenericGA.h"
#include "random.h"
#include "Utils.h"
#include "hps_0.h"
#include <time.h>

// sopc-create-header-files "./testeio.sopcinfo" --single hps_0.h --module hps_0

#define REG_BASE 0xFF200000
#define REG_SPAN 0x00200000

typedef std::array<bool, 16> Func;

//TODO: redefinir os parâmetros da estrutura
#define NUM_IN 2
#define NUM_OUT 1
#define MUTATION_RATE 0.15
#define LAMBDA 4
#define MAX_GENERATIONS 10000

// Circuit parameters
#define CIRCUIT_ROW_COUNT 2
#define CIRCUIT_COLUMN_COUNT 2

using namespace std::chrono; // PARA VERIFICAR TEMPO DE PROCESSAMENTO

std::vector<std::tuple<std::bitset<8>, std::bitset<8>, std::bitset<8>>>
    inputOutputValidSequences() {
    return map([](std::tuple<const char*, const char*, const char*> s)
            { return std::make_tuple(std::bitset<8>(std::get<0>(s)), std::bitset<8>(std::get<1>(s)), std::bitset<8>(std::get<2>(s))); },
            std::vector<std::tuple<const char*, const char*, const char*>>{

				//AND
				std::make_tuple("00000000","00000000", "00000001"),
				std::make_tuple("00000001","00000000", "00000001"),
				std::make_tuple("00000010","00000000", "00000001"),
				std::make_tuple("00000011","00000001", "00000001")
            });
}

//TODO: checar todas as referências do tipo cell.inputs ou cell.function
struct Cell {
	std::array<bool,16> function;
};

struct GeneticParams {
	unsigned int r, c, numIn, numOut, leNumIn;
};

struct Chromosome {
	std::vector<std::vector<Cell>> cells;
	std::vector<unsigned int> outputs;
};

template <typename T, typename U>
std::tuple<T, U> tup(T t, U u) {
    return std::make_tuple(t, u);
}

Chromosome  makeChromosome
		( std::vector<std::vector<Cell>> cells
        , std::vector<unsigned int> outputs
		) {
	Chromosome res;
	res.cells = cells;
	res.outputs = outputs;
	return res;
}

// bool vector has at most 32 elements
uint32_t convertUnit(std::vector<bool> v) {
    uint32_t res = 0;
    for (unsigned int i = 0; i < v.size(); i++) {
        res |= (v[i] << i);
    }
    return res;
}

// bits range = [0..32]
std::vector<bool> convertToBits(uint32_t val, uint32_t bits) {
    std::vector<bool> res;
    for (unsigned int i = 0; i < bits; i++) {
        res.push_back((val & (1 << i)) >> i);
    }
    return res;
}

std::vector<uint32_t> convertToPacked(std::vector<bool> v) {
    return map(convertUnit, chunksOf(v, 32));
}

//TODO: adequar essa funcao de serializacao
std::vector<bool> serializeCell(Cell cell) {
	std::vector<bool> result;

	for(auto f : cell.function){
		result.push_back(f);
	}

	return result;
}

std::vector<bool> rawSerialize(GeneticParams params, Chromosome chrom) {
 	std::vector<uint32_t> result;

	auto numPinos = params.c * params.r;
	auto bitsPinos = (uint32_t) ceil(log2(numPinos));

	std::vector<bool> totalBits;

	for (unsigned int i = 0; i < params.r; i++) {
		for (unsigned int j = 0; j < params.c; j++) {
			auto cellBits = serializeCell(chrom.cells[i][j]);
			totalBits.insert(totalBits.end(), cellBits.begin(), cellBits.end());
		}
	}

	auto outBits = map([=](unsigned int out) {
		return convertToBits(out, bitsPinos);
	}, chrom.outputs);

	for (auto out : outBits) {
        totalBits.insert(totalBits.end(), out.begin(), out.end());
	}

	return totalBits;
}

std::vector<uint32_t> serialize(GeneticParams params, Chromosome chrom) {
	return convertToPacked(rawSerialize(params, chrom));
}

//TODO: Mudar os parametros de toda makeCell pelo codigo
Cell makeCell(std::array<bool,16> func) {
	Cell res;
	res.function = func;
	return res;
}

uint32_t  firstBitOfEachOutput(std::vector<std::bitset<8>> answer) {
	uint32_t result = 0;
	for (unsigned int i = 0; i < answer.size(); i++) {
		result = result + (answer[i][0] << i);
	}
	return result;
}

std::tuple<unsigned int, unsigned int> indexToCoordinate(unsigned int index, unsigned int r) {
	return std::make_tuple(index % r, index / r);
}

unsigned int rawCoordinateToIndex(std::tuple<unsigned int, unsigned int> coordinate, unsigned int r) {
	return std::get<0>(coordinate) + r * std::get<1>(coordinate);
}

unsigned int coordinateToIndex
    ( std::tuple<unsigned int, unsigned int> coordinate
    , unsigned int r
    , unsigned int numIn
    ) {
    return numIn + rawCoordinateToIndex(coordinate, r);
}

unsigned int fitInLargerIndex
    ( unsigned int index
    , unsigned int oldR
    , unsigned int newR
    , unsigned int oldNumIn
    , unsigned int numIn
    , unsigned int lineOffset
    ) {
	if (index >= oldNumIn) {
		index -= oldNumIn;
		auto coord = indexToCoordinate(index, oldR);
		std::get<0>(coord) += lineOffset;
		return coordinateToIndex(coord, newR, numIn);
	}
	return index;
}


//=====================================================================
//TODO: modificar as funcoes que mostram o cromossomo, abaixo
std::string showCell(Cell cell) {

	std::string s = "";

	for(bool f : cell.function){
		f ? s += '1' : s += '0'; 
	}
	
    return s;
}

std::string showInt(unsigned int i) {
	std::ostringstream ss;
	ss << i;
	return ss.str();
}

std::string showInput(GeneticParams params, unsigned int i) {
	std::string s;

	auto col = i % params.c;
	auto row = i / params.c;
	s += "(";
	s += showInt(row);
	s += ", ";
	s += showInt(col);
	s += ")";

	return s;
}

std::vector<std::string> format(std::vector<std::string> vs) {
	auto sizes = map([=](std::string s) { return s.size(); }, vs);
	auto max = std::max_element(sizes.begin(), sizes.end());
	return map([=](std::string s) {
		auto origSize = s.size();
		for (unsigned int i = 0; i < *max - origSize; i++) {
			s += ' ';
		}
		return s;
	}, vs);
}

std::string showChromosome(GeneticParams params, Chromosome chrom) {
	std::string s;
	s += "Cells:\n";
	auto cellsS =
			map([=](std::vector<Cell> row) {
		return map([=](Cell c) { return showCell(c); }, row);
	}, chrom.cells);

	auto tCellsS = transpose(cellsS);

	auto tFormatted = map
			([=](std::vector<std::string> column) {
		return format(column);
	}, tCellsS);

	auto formatted = transpose(tFormatted);

	for (unsigned int i = 0; i < formatted.size(); i++) {
		for (unsigned int j = 0; j < formatted[0].size(); j++) {
			s += formatted[i][j];
			s += "   ";
		}
		s += '\n';
	}

	s += "Outputs:\n";
	for (unsigned int i = 0; i < chrom.outputs.size(); i++) {
		s += showInput(params, chrom.outputs[i]);
		s += ' ';
	}
	s += '\n';
	return s;
}

//===========================================================================

// Side-effectful!! Writes to FPGA ports.
void sendVectorToFPGA(std::vector<uint32_t> vec, void* fpgaMemory) {
	std::vector<int> segmentAddrs = {
		CHROM_SEG_0_BASE, CHROM_SEG_1_BASE, CHROM_SEG_2_BASE, CHROM_SEG_3_BASE, CHROM_SEG_4_BASE,
		CHROM_SEG_5_BASE, CHROM_SEG_6_BASE, CHROM_SEG_7_BASE, CHROM_SEG_8_BASE, CHROM_SEG_9_BASE,
		CHROM_SEG_10_BASE, CHROM_SEG_11_BASE, CHROM_SEG_12_BASE, CHROM_SEG_13_BASE, CHROM_SEG_14_BASE,
		CHROM_SEG_15_BASE, CHROM_SEG_16_BASE, CHROM_SEG_17_BASE, CHROM_SEG_18_BASE, CHROM_SEG_19_BASE,
		CHROM_SEG_20_BASE, CHROM_SEG_21_BASE, CHROM_SEG_22_BASE, CHROM_SEG_23_BASE, CHROM_SEG_24_BASE,
		CHROM_SEG_25_BASE, CHROM_SEG_26_BASE, CHROM_SEG_27_BASE, CHROM_SEG_28_BASE, CHROM_SEG_29_BASE,
		CHROM_SEG_30_BASE
	};

	// Sending serialized chromosome to FPGA.
	for (unsigned int i = 0; i < vec.size(); i++) {
		void* segmentAddr = (uint8_t*) fpgaMemory + segmentAddrs[i];
		*(uint32_t*) segmentAddr = vec[i];
	}
}

void sendChromosomeToFPGA(Chromosome chromosome, GeneticParams params, void* fpgaMemory) {
	sendVectorToFPGA(serialize(params, chromosome), fpgaMemory);
}

uint32_t sendVectorAndGetErrorSum
    ( std::vector<uint32_t> vec
    , void* fpgaMemory
    ) { //COLOCAR PARA MEDIR TEMPO AQUI???

	high_resolution_clock::time_point t5 = high_resolution_clock::now(); //MARCADOR INICIAL DE TEMPO

    std::vector<int> errorSumAddrs = {
        ERROR_SUM_0_BASE, ERROR_SUM_1_BASE, ERROR_SUM_2_BASE, ERROR_SUM_3_BASE,
        ERROR_SUM_4_BASE, ERROR_SUM_5_BASE, ERROR_SUM_6_BASE, ERROR_SUM_7_BASE
    };

    void* doneProcessingFeedbackAddr = (uint8_t*) fpgaMemory + DONE_PROCESSING_FEEDBACK_BASE;
    *(uint32_t*) doneProcessingFeedbackAddr = 0;
	

    void* readyToProcessAddr = (uint8_t*) fpgaMemory + READY_TO_PROCESS_BASE;
    while ((*(uint32_t*) readyToProcessAddr) != 1){
		// Esse print esta aqui porque o otimizador do g++ faz o programa
        // ficar preso num loop infinito se isso nao estiver aqui.
        std::cout << "";
    }

    sendVectorToFPGA(vec, fpgaMemory);

    void* startProcessingAddr = (uint8_t*) fpgaMemory + START_PROCESSING_CHROM_BASE;
    *(uint32_t*) startProcessingAddr = 1;

    void* doneProcessingAddr = (uint8_t*) fpgaMemory + DONE_PROCESSING_CHROM_BASE;
    while ((*(uint32_t*) doneProcessingAddr) != 1) {
        // Esse print esta aqui porque o otimizador do g++ faz o programa
        // ficar preso num loop infinito se isso nao estiver aqui.
        std::cout << "";
    }
		// FIM
		high_resolution_clock::time_point t6 = high_resolution_clock::now(); //MARCADOR FINAL DE TEMPO
		auto tempCrome = duration_cast<microseconds>( t6 - t5 ).count();
		// std::cout << " TempoCromossomo:" << tempCrome << "microssegundos" << std::endl;

    uint32_t chromErrorSum = 0;
    for (auto addr : errorSumAddrs) {
        void* chromErrorSumAddr = (uint8_t*) fpgaMemory + addr;
        chromErrorSum += *(uint32_t*) chromErrorSumAddr;
    }

    // Esse print esta aqui porque o otimizador do g++ faz o programa
    // ficar preso num loop infinito se isso nao estiver aqui.
    std::cout << "";

    *(uint32_t*) startProcessingAddr = 0;
    *(uint32_t*) doneProcessingFeedbackAddr = 1;

    return chromErrorSum;

}

uint32_t sendChromosomeAndGetErrorSum
    ( Chromosome chrom
    , GeneticParams params
    , void* fpgaMemory
    ) {
    auto serial = serialize(params, chrom);
    return sendVectorAndGetErrorSum(serial, fpgaMemory);
}

std::function<double(Chromosome)>
		makeFPGAFitnessFunc
			( GeneticParams params
			, void* fpgaMemory
            ) {
	return [=](Chromosome chrom) {
	    return sendChromosomeAndGetErrorSum(chrom, params, fpgaMemory);
	};
}

std::function<double(Chromosome)>
		makeGrowingFPGAFitnessFunc
			( GeneticParams params
            , GeneticParams largerParams
			, void* fpgaMemory
            ) {
	return [=](Chromosome chrom) {
	    auto largerChrom = chrom;
        return sendChromosomeAndGetErrorSum(largerChrom, largerParams, fpgaMemory);
	};
}

std::function<double(Chromosome)>
		makeFPGASeparateFitnessFunc
			( GeneticParams params
			, unsigned int outputNum
			, unsigned int numOutputs
			, void* fpgaMemory
            ) {
	return [=](Chromosome chrom) {
		auto largerParams = params;
		largerParams.r = params.r * numOutputs;
		largerParams.numOut = numOutputs;

		return sendChromosomeAndGetErrorSum(chrom, largerParams, fpgaMemory);
	};
}

std::array<bool,16> arrayFromInt(int r){
	std::array<bool,16> a;

	a[0] = r & 0x0001;
	a[1] = (r & 0x0002) >> 1;
	a[2] = (r & 0x0004) >> 2;
	a[3] = (r & 0x0008) >> 3;
	a[4] = (r & 0x0010) >> 4;
	a[5] = (r & 0x0020) >> 5;
	a[6] = (r & 0x0040) >> 6;
	a[7] = (r & 0x0080) >> 7;
	a[8] = (r & 0x0100) >> 8;
	a[9] = (r & 0x0200) >> 9;
	a[10] = (r & 0x0400) >> 10;
	a[11] = (r & 0x0800) >> 11;
	a[12] = (r & 0x1000) >> 12;
	a[13] = (r & 0x2000) >> 13;
	a[14] = (r & 0x4000) >> 14;
	a[15] = (r & 0x8000) >> 15;

	return a;
}

RNGFUNC(Func) randomFunc() {
	return rmap<random_type, Func>([](random_type r) {
		return arrayFromInt(r);
	}, getRandom());
}

RNGFUNC(unsigned int)  randomOutput(GeneticParams params, unsigned int c) {
	return bind
			( getRandom()
            , [=](random_type rand) {
		return pure(rand % (params.numIn + params.r * params.c));
	});
}

RNGFUNC(std::vector<unsigned int>)  mutateOutput
		( std::vector<unsigned int> outputs
        , unsigned int pointToMutate
        , GeneticParams params) {
	return bind
			( randomOutput(params, params.c)
            , [=](unsigned int newOut) mutable {
		outputs[pointToMutate] = newOut;
		return pure(outputs);
	});
}

//TODO: Gerar o novo cromossomo aleatorio e ajustar a mutacao
RNGFUNC(Cell)  randomCell()
{
	return bind(randomFunc(), [=](Func randFunc) {
		return pure(makeCell(randFunc));
	});
}

RNGFUNC(Cell) mutateCell(Cell cell){
	return bind(getRandom(), [=](random_type rand) {
		Cell cellMut = cell;
		cellMut.function[rand % 15] = !cellMut.function[rand % 15];
		return pure(cellMut);
	});
}

RNGFUNC(std::vector<std::vector<Cell>>)  mutateGrid
		( std::vector<std::vector<Cell>> grid
        , random_type pointToMutate
		, GeneticParams params
		) {
	return bind(mutateCell(grid[pointToMutate/params.c][pointToMutate%params.c]), [=](Cell cell) {
		auto newGrid = grid;
		newGrid[pointToMutate/params.c][pointToMutate%params.c] = cell;
		return pure(newGrid);
    });
}

std::function<RNGFUNC(Chromosome)(Chromosome)>
	makeMutation(GeneticParams params, float mutationPercent) {


    auto totalElements = params.r * params.c * 16 + params.numOut;
    auto elementsToMutate = std::ceil(totalElements * mutationPercent);

	return [=](Chromosome chrom) {
		return bind(sequence(replicate(elementsToMutate, getRandom())),
				[=](std::vector<random_type> rands) mutable {
			return foldM([=](Chromosome c, random_type rand) mutable {
                auto pointToMutate = rand % (params.c * params.r + params.numOut);
                if (pointToMutate < params.c * params.r) {
                    return bind
                            ( mutateGrid(chrom.cells, pointToMutate, params)
                            , [=](std::vector<std::vector<Cell>> newGrid) mutable {
                        c.cells = newGrid;
                        return pure(c);
                    });
                } else {
                    return bind
                            ( mutateOutput(chrom.outputs, pointToMutate - (params.c * params.r), params)
                            , [=](std::vector<unsigned int> newOuts) mutable {
                        c.outputs = newOuts;
                        return pure(c);
                    });
                }
			}, chrom, rands);
		});
	};

}

RNGFUNC(std::vector<Cell>)  randomColumn(GeneticParams params) {
	return sequence(replicate(params.r, randomCell()));
}

RNGFUNC(std::vector<std::vector<Cell>>)  randomCells(GeneticParams params) {
	return bind(mapM([=](unsigned int c) {
		return randomColumn(params);
	}, vectorFromTo(0, params.c)), [=](std::vector<std::vector<Cell>> grid) {
		return pure(transpose(grid));
	});
}

RNGFUNC(std::vector<unsigned int>)  randomOutputs(GeneticParams params) {
	return sequence(replicate(params.numOut, randomOutput(params, params.c)));
}

RNGFUNC(Chromosome)  randomChrom(GeneticParams params) {
	return bind(randomOutputs(params), [=](std::vector<unsigned int> randomOuts) {
		return bind(randomCells(params), [=](std::vector<std::vector<Cell>> randomCs) {
			return pure(makeChromosome(randomCs, randomOuts));
		});
	});
}

void dumpMemoryToFile(void* fpgaMem, uint32_t base, std::string fileName) {
    std::ofstream out;
    out.open(fileName.c_str());

    uint8_t* memBase = (uint8_t*) fpgaMem + base;
    for (unsigned int i = 0; i < 32768; i++) {
        auto val = *(uint32_t*) (memBase + (i * 4));
        out << std::hex << std::setfill('0') << std::setw(8) << val << std::endl;
    }

    out.close();
}

template <typename F>
std::function<bool(GAState<Evaluated<Chromosome>>)>
    makeCorrectTermination(GeneticParams currentParams, GeneticParams finalParams, F fitness, void* fpgaMem) {

    static_assert(std::is_convertible<F, std::function<double(Chromosome)>> ::value,
                "any's function must be of type Chromosome -> double");

    return [=](GAState<Evaluated<Chromosome>> state) {
        if (state.population[0].score == 0) {
            auto verifiedScore = fitness(state.population[0].value);
            if (verifiedScore == 0) {
                printf("%d %g\n", state.generation, state.population[0].score);
                return false;
            }
            auto fileName = std::string("Failed_") + showInt(currentParams.r) + "_" + showInt(state.generation) + ".txt";
            std::cout << "Failed re-verification. Writing to file " << fileName << std::endl;

            std::ofstream out;
            out.open(fileName.c_str());

            out << showChromosome(currentParams, state.population[0].value) << std::endl;
            auto s = map([](bool b) { return b ? '1' : '0'; }, rawSerialize(currentParams, state.population[0].value));
            std::reverse(s.begin(), s.end());
            out << "Normal bitstring:" << std::endl;
            out << std::string(s.begin(), s.end()) << std::endl << std::endl;
	
            auto largerChrom = state.population[0].value;
            s = map([](bool b) { return b ? '1' : '0'; }, rawSerialize(finalParams, largerChrom));
            std::reverse(s.begin(), s.end());
            out << "Larger bitstring:" << std::endl;
            out << std::string(s.begin(), s.end()) << std::endl << std::endl;
            out.close();

            dumpMemoryToFile(fpgaMem, TWO_PORT_MEM_BASE, std::string("Dump_") + showInt(currentParams.r) + "_" + showInt(state.generation) + ".txt");
            dumpMemoryToFile(fpgaMem, TWO_PORT_MEM_CORRECT_BASE, std::string("Correct_Dump_") + showInt(currentParams.r) + "_" + showInt(state.generation) + ".txt");

            return true;
        }
        if (state.generation % 1 == 0) {	//MOSTRAR O NUMERO DE CADA GERACAO
            printf("%d %g\n", state.generation, state.population[0].score);
        }
        return state.generation < MAX_GENERATIONS;
    };
}

bool correctTermination(GAState<Evaluated<Chromosome>> state) {
    if (state.generation % 10 == 0 || state.population[0].score == 0) {
        printf("%d %g\n", state.generation, state.population[0].score);
    }
	return state.generation < 50000 && state.population[0].score > 0;
}

bool  optimizingTermination(GAState<Evaluated<Chromosome>> state) {
	printf("%d %g\n", state.generation, state.population[0].score);
	return state.generation < 50000;
}

template <typename F>
RNGFUNC(GAState<Evaluated<Chromosome>>)
	GARoutineWithInitial
		( GeneticParams params
		, GeneticParams finalParams
		, F fitnessFunc
		, Chromosome initial
		, void* fpgaMem
        ) {

    static_assert(std::is_convertible<F, std::function<double(Chromosome)>> ::value,
                "fpgaGARoutine's fitness function must be of type T -> bool");

    auto mutationFunc = makeMutation(params, MUTATION_RATE);
    auto strategy = lambdaPlusN<Chromosome>(fitnessFunc, mutationFunc, LAMBDA); //MEDIA 1+LAMBDA AQUI??
    auto gaFunc = makeGAFunction<Evaluated<Chromosome>>(strategy);
    auto termination = makeCorrectTermination(params, finalParams, fitnessFunc, fpgaMem);

    GAState<Evaluated<Chromosome>> init;
    init.generation = 0;
    init.population = { makeEvaluated(initial, fitnessFunc(initial)) };

    return iterateWhileM(termination, gaFunc, init);
}

RNGFUNC(std::vector<GAState<Evaluated<Chromosome>>>)
	fpgaSeparateGARoutine
		( GeneticParams params
		, unsigned int totalNumOutputs
        , void* fpgaMemory) {
	return mapM([=](unsigned int outputNum) {
        return bind
        		( randomChrom(params)
                , [=](Chromosome randomChromosome) {
            return GARoutineWithInitial
                       ( params
                       , params
                       , makeFPGASeparateFitnessFunc(params, outputNum, totalNumOutputs, fpgaMemory)
                       , randomChromosome
                       , fpgaMemory
                       );
        });
	}, vectorFromTo(0, totalNumOutputs));
}

RNGFUNC(GAState<Evaluated<Chromosome>>)
	fpgaGARoutine
		( GeneticParams params
        , void* fpgaMemory
        ) {
    return bind
            ( randomChrom(params)
            , [=](Chromosome randomChromosome) {
        return GARoutineWithInitial
                   ( params
                   , params
                   , makeFPGAFitnessFunc(params, fpgaMemory)
                   , randomChromosome
                   , fpgaMemory
                   );
    });
}

RNGFUNC(GAState<Evaluated<Chromosome>>)
	fpgaGrowingGARoutine
		( GeneticParams params
        , GeneticParams largerParams
        , void* fpgaMemory
        ) {
    return bind
            ( randomChrom(params)
            , [=](Chromosome randomChromosome) {
        return GARoutineWithInitial
                   ( params
                   , largerParams
                   , makeGrowingFPGAFitnessFunc(params, largerParams, fpgaMemory)
                   , randomChromosome
                   , fpgaMemory
                   );
    });
}

void* openFPGAMemory() {
	int fd = open("/dev/mem", O_RDWR | O_SYNC);
	if (!fd) {
		std::cout << "Could not open /dev/mem." << std::endl;
		exit(1);
	}

	void* virtualBase = mmap(NULL, REG_SPAN, PROT_READ | PROT_WRITE, MAP_SHARED, fd, REG_BASE);
	if (!virtualBase) {
		std::cout << "Could not map physical memory into virtual." << std::endl;
		exit(1);
	}
	return virtualBase;
}

std::vector<GeneticParams> growingGenParamVec(GeneticParams baseParams, unsigned int initialR, unsigned int finalR) {
    return map([=](unsigned int cur) {
        GeneticParams p = baseParams;
        p.r = cur;
        return p;
    }, vectorFromTo(initialR, finalR + 1));
}

int main() {
	GeneticParams params;
	params.r = 1; // Para cada solucao individual.
	params.c = CIRCUIT_COLUMN_COUNT;
	params.numIn = NUM_IN;
	params.numOut = NUM_OUT;
	params.leNumIn = 2;

	srand(time(NULL));
	auto initialRng = rand();

	printf("Seed: %d\n", initialRng);

	auto fpgaMem = openFPGAMemory();

	auto io = inputOutputValidSequences();

	// Send the size of the sequences to be processed
	void* seqToProcAddr = (uint8_t*) fpgaMem + SEQUENCES_TO_PROCESS_BASE;
	*(uint32_t*) seqToProcAddr = io.size();

	//Seleciona a pio 
	void* inputAddress = (uint8_t*) fpgaMem + INPUT_SEQUENCE_0_BASE;
	void* outputAddress = (uint8_t*) fpgaMem + EXPECTED_OUTPUT_0_BASE;
	void* validOutputAddress = (uint8_t*) fpgaMem + VALID_OUTPUT_0_BASE;

	//Sinais de controle serial
	//Entrada Processador C++
	void* nextSampleAddress = (uint8_t*) fpgaMem + NEXTSAMPLE_BASE;

	//Saida processador C++
	void* sampleIndexAddress = (uint8_t*) fpgaMem + SAMPLEINDEX_BASE;
	void* preparingNextSampleAddress = (uint8_t*) fpgaMem + PREPARINGNEXTSAMPLE_BASE;
	void* writeSampleAddress = (uint8_t*) fpgaMem + WRITESAMPLE_BASE;

	//Inicializando saidas
	*(uint32_t*) preparingNextSampleAddress = 0;
	*(uint32_t*) sampleIndexAddress = 0;		//entrada na maquina de estados
	*(uint32_t*) writeSampleAddress = 0;	

	high_resolution_clock::time_point t1 = high_resolution_clock::now(); //INICIO DO MARCADOR DE TEMPO PARA ENVIO DOS INDIVIDUOS A FPGA

	for (unsigned int i = 0; i < io.size(); i += 4) {
	    for(unsigned int j = 0; j < 4; j++){
			auto seg = j % 4; //precisa pra dividir os 32 bits em 4 de 8

			unsigned int linha = i+j;
			
			printf("amostra %d\n", i+j);
			*(uint32_t*) preparingNextSampleAddress = 1; //Manda FPGA sair do IDLE
			while((*(uint32_t*) nextSampleAddress) == 1){ //Espera FPGA sair do IDLE
				// Esse print esta aqui porque o otimizador do g++ faz o programa
				// ficar preso num loop infinito se isso nao estiver aqui.
				std::cout << "";
			}

			*(uint32_t*) preparingNextSampleAddress = 0;

			//Transforma o bitset da seq em 32 bits e faz shift
			uint32_t inputVal = std::get<0>(io[linha]).to_ulong() << (8 * seg);
			uint32_t outputVal = std::get<1>(io[linha]).to_ulong() << (8 * seg);
			uint32_t validOutputVal = std::get<2>(io[linha]).to_ulong() << (8 * seg);
			
			//Insere a seq no pio diretamente (1 vez) ou com or(|)
			*(uint32_t*) inputAddress = seg == 0 ? inputVal : *(uint32_t*) inputAddress | inputVal;
			*(uint32_t*) outputAddress = seg == 0 ? outputVal : *(uint32_t*) outputAddress | outputVal;
			*(uint32_t*) validOutputAddress = seg == 0 ? validOutputVal : *(uint32_t*) validOutputAddress | validOutputVal;		
		}

		//Escreve o indice
		*(uint32_t*) sampleIndexAddress = i/4;

		printf("ESCRITA INDICE %d\n",i/4);
		std::cout << "input : " << std::bitset<32>(*(uint32_t*) inputAddress) << "\n";
		std::cout << "output: " << std::bitset<32>(*(uint32_t*) outputAddress) << "\n";
		std::cout << "valid : " << std::bitset<32>(*(uint32_t*) validOutputAddress) << "\n";
 
		//Manda a FPGA escrever
		*(uint32_t*) writeSampleAddress = 1;


		//Aguarda fim da escrita e FPGA retornar para IDLE
		while(*(uint32_t*) nextSampleAddress == 0){
			// Esse print esta aqui porque o otimizador do g++ faz o programa
			// ficar preso num loop infinito se isso nao estiver aqui.
			std::cout << "";
		}
		
		*(uint32_t*) writeSampleAddress = 0;

	}

	 //FINAL DO MARCADOR DE TEMPO DO ENVIO DE INDIVIDUOS A FPGA
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	auto tempGasto = duration_cast<microseconds>( t2 - t1 ).count();

	// std::cout << " Tempo de envio para FPGA:" << tempGasto << "microssegundos" << std::endl;
		printf("Inicio da evolucao!!!\n\n\n");

	/* Send a raw chromosome and evaluate once

	std::string bitstring = "000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000101001000000000001001001001001";
	std::reverse(bitstring.begin(), bitstring.end());
	auto data = std::vector<char>(bitstring.begin(), bitstring.end());
	auto fun = [](char c) {
	    return c == '1' ? true : false;
	};
	std::cout << sendVectorAndGetErrorSum(convertToPacked(map(fun, data)), fpgaMem) << std::endl;

	return 0;
	 */

	GeneticParams finalParams = params;
	finalParams.r = CIRCUIT_ROW_COUNT;
	finalParams.numIn = NUM_IN;
	finalParams.numOut = NUM_OUT;

	auto b2c = [](bool b) {
	        return b ? '1' : '0';
	    };

	auto solution = whileFoldM([=](GeneticParams currentParams, GAState<Evaluated<Chromosome>> finalSolution) {
	    if (finalSolution.population[0].score > 0) {
	        std::cout << "Increasing number of rows by one, new value: " << currentParams.r + 1 << std::endl;
	        return true;
	    }

        printf("Solution:\n");
        printf("%s\n", showChromosome(currentParams, finalSolution.population[0].value).c_str());
        auto s = map(b2c, rawSerialize(currentParams, finalSolution.population[0].value));
        std::reverse(s.begin(), s.end());
        std::cout << "Normal bitstring:" << std::endl;
        std::cout << std::string(s.begin(), s.end()) << std::endl << std::endl;

        auto largerChrom = finalSolution.population[0].value;
        s = map(b2c, rawSerialize(finalParams, largerChrom));
        std::reverse(s.begin(), s.end());
        std::cout << "Larger bitstring:" << std::endl;
        std::cout << std::string(s.begin(), s.end()) << std::endl << std::endl;

        auto fit = makeGrowingFPGAFitnessFunc(currentParams, finalParams, fpgaMem)(finalSolution.population[0].value);
        printf("Fitness recalculated: %g\n", fit);

        dumpMemoryToFile(fpgaMem, TWO_PORT_MEM_BASE, "Dump_Recalculated.txt");

	    return false;
	},
	[fpgaMem, finalParams](GeneticParams params) {
	    return fpgaGrowingGARoutine(params, finalParams, fpgaMem);
	}, growingGenParamVec(params, CIRCUIT_ROW_COUNT, CIRCUIT_ROW_COUNT));

	evalState(solution, initialRng);

	return 0;
}
