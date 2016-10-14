/*
 * TSPsolver.h
 *
 *  Created on: 20/09/2016
 *      Author: romanelli
 */

#ifndef TSPSOLVER_H_
#define TSPSOLVER_H_

#include <string>

namespace TSPsolverOpcao {
	const int OpcaoBuscaLocalConvencional = 0;
	const int OpcaoBuscaLocalRapidaArestasAleatorias = 1;
	const int OpcaoBuscaLocalRapidaArestasMenoresPrimeiro = 2;
	const int OpcaoBuscaLocalRapidaArestasMaioresPrimeiro = 3;
	const std::string StrOpcaoBuscaLocalConvencional = "BLC";
	const std::string StrOpcaoBuscaLocalRapidaArestasAleatorias = "BLRAle";
	const std::string StrOpcaoBuscaLocalRapidaArestasMenoresPrimeiro = "BLRMen";
	const std::string StrOpcaoBuscaLocalRapidaArestasMaioresPrimeiro = "BLRMai";

	const int OpcaoPrimeiroAprimorante = 0;
	const int OpcaoMelhorAprimorante = 1;
}

class TSPsolver {
private:
	int numVertices;
	double** pesos;
	int** penalidade;
	double lambda;
	int iteracoes;
	bool* featureSubNeighborhoodActivation;
	int numSubNeighborhoodsActive;
	int** sucessor;
	int opcao;
public:
	TSPsolver(int nv, double** c, double lambda, int iteracoes, int opcao);
	double funcaoCustoSolucao(int* solucao);
	double funcaoCustoSolucaoAumentada(int* solucao);
	double expressaoUtilidade(int* solucao, int indCaracteristica);
	int* resolver();
	int* guidedLocalSearch();
	int* guidedLocalSearchWithFastLocalSearch();
	int* localSearch(int* solucaoAtual, int* melhorSolucao,
			bool usarFuncaoCustoAumentada);
	int* fastLocalSearch(int* solucaoAtual, int* melhorSolucao, int opcaoAprimorante,
			bool usarFuncaoCustoAumentada);
	int* gerarSolucaoArbitraria();
	double obterPesoAresta(int origem, int destino);
	int obterPenalidadeAresta(int origem, int destino);
	void calcularMatrizSucessoresOrdenada(bool ordemCrescente);
	void calcularMatrizSucessoresAleatoria();
	virtual ~TSPsolver();
};

#endif /* TSPSOLVER_H_ */
















