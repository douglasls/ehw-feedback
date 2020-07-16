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
#include <stdint.h>
#include <stdio.h>
#include <cassert>
#include "GenericGA.h"
#include "random.h"
#include "Utils.h"
#include "hps_0.h"


std::vector<std::tuple<std::bitset<8>, std::bitset<8>, std::bitset<8>>>
inputOutputValidSequences() {
return map([](std::tuple<const char*, const char*, const char*> s)
        { return std::make_tuple(std::bitset<8>(std::get<0>(s)), std::bitset<8>(std::get<1>(s)), std::bitset<8>(std::get<2>(s))); },
        std::vector<std::tuple<const char*, const char*, const char*>>{

            
            //SOMADOR COMPLETO 1BIT 3i_2o
            std::make_tuple("00000000", "00000000", "00000011"),
            std::make_tuple("00000001", "00000000", "00000011"),
            std::make_tuple("00000101", "00000010", "00000011"),
            std::make_tuple("00000100", "00000001", "00000011"),

            std::make_tuple("00000000", "00000000", "00000011"),
            std::make_tuple("00000100", "00000001", "00000011"),
            std::make_tuple("00000110", "00000010", "00000011"),
            std::make_tuple("00000010", "00000001", "00000011"),

            std::make_tuple("00000011", "00000010", "00000011"),
            std::make_tuple("00000010", "00000001", "00000011"),
            std::make_tuple("00000000", "00000000", "00000011"),
            std::make_tuple("00000001", "00000001", "00000011"),

            std::make_tuple("00000011", "00000010", "00000011"),
            std::make_tuple("00000111", "00000011", "00000011"),
            std::make_tuple("00000011", "00000010", "00000011"),
            std::make_tuple("00000001", "00000001", "00000011"),

            std::make_tuple("00000101", "00000010", "00000011"),
            std::make_tuple("00000111", "00000011", "00000011"),
            std::make_tuple("00000101", "00000010", "00000011"),
            std::make_tuple("00000100", "00000001", "00000011"),

            std::make_tuple("00000110", "00000010", "00000011"),
            std::make_tuple("00000111", "00000011", "00000011"),
            std::make_tuple("00000110", "00000010", "00000011"),
            std::make_tuple("00000010", "00000001", "00000011"),

            std::make_tuple("00000000", "00000000", "00000011"),
            std::make_tuple("00000000", "00000000", "00000000"), //FIO
            std::make_tuple("00000000", "00000000", "00000000"), //FIO
            std::make_tuple("00000000", "00000000", "00000000")  //FIO


        });
}

int main(){
    auto io = inputOutputValidSequences();

    std::cout << typeid(std::get<0>(io[0])).name() << std::endl;

    //get<coluna>( io[linha] )
    std::cout << std::get<0>(io[10]) << std::endl;

    return 0;
}