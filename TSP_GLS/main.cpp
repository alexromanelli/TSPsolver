/*
 * main.cpp
 *
 *  Created on: 20/09/2016
 *      Author: romanelli
 */

#include <iostream>
#include <cstdio>
#include <fstream>
#include <string>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <algorithm>
#include <functional>
#include <cctype>
#include <locale>
#include <cmath>
#include <time.h>
#include <sys/time.h>

#include "TSPsolver.h"

double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}

// trim from start
static inline std::string &ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
            std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(),
            std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
}

// trim from both ends
static inline std::string &trim(std::string &s) {
    return ltrim(rtrim(s));
}

void split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (getline(ss, item, delim)) {
    	if (item != "")
    		elems.push_back(item);
    }
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

void lerMatrizDePesos(std::ifstream* arqInstancia, int numCidades, double** pesos) {
	std::string linha;
	for (int i = 0; i < numCidades - 1; i++) {
		getline(*arqInstancia, linha);
		std::vector<std::string> spesos = split(trim(linha), ' ');
		for (unsigned int j = 0; j < spesos.size(); j++) {
			double peso = atof(spesos[j].c_str());
			pesos[i][i + j + 1] = peso;
		}
	}
}

double calcularDistancia(double x1, double y1, double x2, double y2) {
	double dist = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
	return dist;
}

void lerCoordenadas(std::ifstream* arqInstancia, int numCidades, double** pesos) {
	std::string linha;
	double** coord = new double*[numCidades];
	for (int i = 0; i < numCidades; i++)
		coord[i] = new double[2]; // 2 : lat, lon
	for (int i = 0; i < numCidades; i++) {
		getline(*arqInstancia, linha);
		std::vector<std::string> snums = split(trim(linha), ' ');
		coord[i][0] = atof(snums[1].c_str());
		coord[i][1] = atof(snums[2].c_str());
	}
	// calcular os pesos (distâncias entre as cidades)
	for (int i = 0; i < numCidades - 1; i++) {
		for (int j = i + 1; j < numCidades; j++) {
			pesos[i][j] = calcularDistancia(coord[i][0], coord[i][1],
					coord[j][0], coord[j][1]);
		}
	}
}

void imprimirPesos(int numCidades, double** pesos) {
	for (int i = 0; i < numCidades; i++) {
		for (int j = 0; j < numCidades; j++) {
			std::printf(" %6.1f ", pesos[i][j]);
		}
		std::printf("\n");
	}
}

int main(int argc, char* argv[]) {
	if (argc == 5) {
		double lambda = atof(argv[1]);
		char* nomeArqInstancia = argv[2];
		std::printf("%s\n", nomeArqInstancia);
		std::string nomeMetodoBL = std::string(argv[3]);

		int opcao = -1;
		if (nomeMetodoBL == TSPsolverOpcao::StrOpcaoBuscaLocalConvencional)
			opcao = TSPsolverOpcao::OpcaoBuscaLocalConvencional;
		else if (nomeMetodoBL == TSPsolverOpcao::StrOpcaoBuscaLocalRapidaArestasAleatorias)
			opcao = TSPsolverOpcao::OpcaoBuscaLocalRapidaArestasAleatorias;
		else if (nomeMetodoBL == TSPsolverOpcao::StrOpcaoBuscaLocalRapidaArestasMaioresPrimeiro)
			opcao = TSPsolverOpcao::OpcaoBuscaLocalRapidaArestasMaioresPrimeiro;
		else if (nomeMetodoBL == TSPsolverOpcao::StrOpcaoBuscaLocalRapidaArestasMenoresPrimeiro)
			opcao = TSPsolverOpcao::OpcaoBuscaLocalRapidaArestasMenoresPrimeiro;

		std::ifstream arqInstancia(nomeArqInstancia);
		int numIteracoes = atoi(argv[4]);

		if (arqInstancia.is_open()) {
			std::string linha = "";
			int numCidades = 0;
			int estado = 0;
			double** pesos;
			while (getline(arqInstancia, linha)) {
				switch (estado) {
				case 0: // procurando número de cidades (dimensão)
					if (linha.substr(0, 9) == "DIMENSION") {
						std::string sDim = linha.substr(11);
						numCidades = atoi(sDim.c_str());
						pesos = new double*[numCidades];
						for (int i = 0; i < numCidades; i++)
							pesos[i] = new double[numCidades];
						estado = 1;
					}
					break;
				case 1: // procurando seção de pesos ou de coordenadas
					if (linha.substr(0, 19) == "EDGE_WEIGHT_SECTION") {
						lerMatrizDePesos(&arqInstancia, numCidades, pesos);
					} else if (linha.substr(0, 18) == "NODE_COORD_SECTION") {
						lerCoordenadas(&arqInstancia, numCidades, pesos);
					}
					break;
				}
			}
			imprimirPesos(numCidades, pesos);
			arqInstancia.close();

			TSPsolver* tspSolver = new TSPsolver(numCidades, pesos, lambda, numIteracoes, opcao);
			int* rota = tspSolver->resolver();

			std::printf("\nMelhor rota encontrada:\n");
			for (int i = 0; i < numCidades; i++) {
				std::printf("%s %d ", (i > 0 ? "," : ""), rota[i] + 1);
			}
			std::printf("\n");
			std::printf("Custo da melhor rota encontrada: %.1f\n", tspSolver->funcaoCustoSolucao(rota));
		} else {
			printf("Erro ao abrir arquivo de instância.\n");
		}
	} else {
		printf("Este programa requer dois parâmetros:\n - o arquivo de instância;\n - o método de busca local, greedy ou fls (fast local search).\n");
	}

	std::cout << "Tempo de execução: " << get_cpu_time() << std::endl;

	std::cout << "Fim da execução." << std::endl;
}


