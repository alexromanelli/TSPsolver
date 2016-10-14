/*
 * TSPsolver.cpp
 *
 *  Created on: 20/09/2016
 *      Author: romanelli
 */

#include "TSPsolver.h"

#include <stdlib.h>
#include <cstdio>
#include <algorithm>
#include <cfloat>
#include <time.h>
#include <string>
#include <cmath>

TSPsolver::TSPsolver(int nv, double** c, double lambda, int iteracoes, int opcao) {
	this->numVertices = nv;
	this->lambda = lambda;
	this->iteracoes = iteracoes;
	this->opcao = opcao;
	this->pesos = new double*[this->numVertices - 1];
	for (int i = 0; i < this->numVertices - 1; i++) {
		this->pesos[i] = new double[this->numVertices];
		for (int j = 0; j < this->numVertices; j++)
			this->pesos[i][j] = c[i][j];
	}
	this->penalidade = new int*[this->numVertices - 1];
	for (int i = 0; i < this->numVertices - 1; i++) {
		this->penalidade[i] = new int[this->numVertices];
	}
	this->featureSubNeighborhoodActivation = new bool[this->numVertices];
}

TSPsolver::~TSPsolver() {
	for (int i = 0; i < this->numVertices - 1; i++) {
		delete[] this->pesos[i];
		delete[] this->penalidade[i];
	}
}

int* TSPsolver::resolver() {
	if (this->lambda == -1) {
		// calcular lambda
		int* solucaoAleatoria = this->gerarSolucaoArbitraria();
		int* s = new int[this->numVertices];
		int* solucao = this->localSearch(solucaoAleatoria, s, false);
		double custo = this->funcaoCustoSolucao(solucaoAleatoria);
		lambda = 0.3 * custo / this->numVertices;
	}
	if (this->opcao == TSPsolverOpcao::OpcaoBuscaLocalConvencional)
		return guidedLocalSearch();
	else
		return guidedLocalSearchWithFastLocalSearch();
}

void ajustaIndices(int* atual, int* sucessor) {
	if (*atual > *sucessor) {
		int t = *atual;
		*atual = *sucessor;
		*sucessor = t;
	}
}

double TSPsolver::funcaoCustoSolucao(int* solucao) {
	double soma = 0;
	for (int i = 0; i < this->numVertices; i++) {
		int atual = solucao[i];
		int sucessor = solucao[(i + 1) % this->numVertices];
		ajustaIndices(&atual, &sucessor);
		double peso = this->pesos[atual][sucessor];
		soma += peso;
	}
	return soma;
}

double TSPsolver::funcaoCustoSolucaoAumentada(int* solucao) {
	double somaCusto = 0;
	double somaRegularizacao = 0;
	for (int i = 0; i < this->numVertices; i++) {
		int atual = solucao[i];
		int sucessor = solucao[(i + 1) % this->numVertices];
		ajustaIndices(&atual, &sucessor);
		double peso = this->pesos[atual][sucessor];
		double penalidade = this->penalidade[atual][sucessor];
		somaCusto += peso;
		somaRegularizacao += penalidade;
	}
	return somaCusto + this->lambda * somaRegularizacao;
}

double TSPsolver::expressaoUtilidade(int* solucao, int indCaracteristica) {
	int atual = solucao[indCaracteristica];
	int sucessor = solucao[(indCaracteristica + 1) % this->numVertices];
	ajustaIndices(&atual, &sucessor);
	double util = this->pesos[atual][sucessor] / (1 + this->penalidade[atual][sucessor]);
	return util;
}

void copiarVetor(int* origem, int* destino, int numElementos) {
	for (int i = 0; i < numElementos; i++) {
		destino[i] = origem[i];
	}
}

int* TSPsolver::guidedLocalSearch() {
	int k = 0;
	int* solucaoInicial = this->gerarSolucaoArbitraria();
	int* melhorSolucao = new int[this->numVertices];
	copiarVetor(solucaoInicial, melhorSolucao, this->numVertices);
	for (int i = 0; i < this->numVertices - 1; i++)
		for (int j = 0; j < this->numVertices; j++)
			this->penalidade[i][j] = 0;

	int* solucaoAtual = solucaoInicial;
	double* util = new double[this->numVertices];

	int* solucaoSucessora;

	while (k < this->iteracoes) {
		std::printf("Iteração: %d\n", k);

		solucaoSucessora = this->localSearch(solucaoAtual, melhorSolucao, true);

		double maxUtil = 0;
		for (int i = 0; i < this->numVertices; i++) {
			util[i] = this->expressaoUtilidade(solucaoSucessora, i);
			if (util[i] > maxUtil)
				maxUtil = util[i];
		}
		for (int i = 0; i < this->numVertices; i++) {
			if (util[i] == maxUtil) {
				int atual = solucaoSucessora[i];
				int sucessor = solucaoSucessora[(i + 1) % this->numVertices];
				ajustaIndices(&atual, &sucessor);
				this->penalidade[atual][sucessor]++;

				std::printf("--> aresta (%d, %d) penalizada: %d\n", atual + 1, sucessor + 1, this->penalidade[atual][sucessor]);
			}
		}

		k++;
	}

	delete[] solucaoAtual;
	delete[] solucaoSucessora;

	return melhorSolucao;
}

void TSPsolver::calcularMatrizSucessoresOrdenada(bool ordemCrescente) {
	for (int i = 0; i < this->numVertices; i++) {
		// iniciar adjacentes com sequência
		for (int j = 0; j < this->numVertices; j++)
			this->sucessor[i][j] = j;

		// ordenar adjacentes por menor distância (observar peso)
		for (int j = 0; j < this->numVertices - 1; j++) {
			double pesoMelhor = i == this->sucessor[i][j] ? DBL_MAX :
					obterPesoAresta(i, this->sucessor[i][j]) + this->lambda * obterPenalidadeAresta(i, this->sucessor[i][j]);
			int indMelhor = j;
			for (int k = j + 1; k < this->numVertices; k++) {
				double pesoK = i == this->sucessor[i][k] ? DBL_MAX :
						obterPesoAresta(i, this->sucessor[i][k]) + this->lambda * obterPenalidadeAresta(i, this->sucessor[i][k]);
				if (ordemCrescente) {
					if (pesoK < pesoMelhor) {
						pesoMelhor = pesoK;
						indMelhor = k;
					}
				} else {
					if (pesoK > pesoMelhor) {
						pesoMelhor = pesoK;
						indMelhor = k;
					}
				}
			}
			// trocar j com indMelhor
			int temp = this->sucessor[i][j];
			this->sucessor[i][j] = this->sucessor[i][indMelhor];
			this->sucessor[i][indMelhor] = temp;
		}
	}
}

void TSPsolver::calcularMatrizSucessoresAleatoria() {
	srand(time(NULL));
	for (int i = 0; i < this->numVertices; i++) {
		// iniciar adjacentes com sequência
		for (int j = 0; j < this->numVertices; j++)
			this->sucessor[i][j] = j;

		// efetuar numVertices trocas aleatórias
		for (int j = 0; j < this->numVertices; j++) {
			int ind1 = rand() % this->numVertices;
			int ind2 = rand() % this->numVertices;
			if (ind2 == ind1)
				ind2 = (ind1 * 3) % this->numVertices;
			// trocar valores das posições ind1 e ind2
			int temp = this->sucessor[i][ind1];
			this->sucessor[i][ind1] = this->sucessor[i][ind2];
			this->sucessor[i][ind2] = temp;
		}
	}
}

int* TSPsolver::guidedLocalSearchWithFastLocalSearch() {
	// criar matriz de sucessores ordenados por distância a partir de um vértice
	this->sucessor = new int*[this->numVertices];
	for (int i = 0; i < this->numVertices; i++) {
		this->sucessor[i] = new int[this->numVertices];
	}

	for (int i = 0; i < this->numVertices; i++) { // for every city
		featureSubNeighborhoodActivation[i] = true; // activate feature neighborhood
	}
	this->numSubNeighborhoodsActive = this->numVertices;

	int k = 0;
	int* solucaoInicial = this->gerarSolucaoArbitraria();
	int* melhorSolucao = new int[this->numVertices];
	copiarVetor(solucaoInicial, melhorSolucao, this->numVertices);
	for (int i = 0; i < this->numVertices - 1; i++)
		for (int j = 0; j < this->numVertices; j++)
			this->penalidade[i][j] = 0;

	int* solucaoAtual = solucaoInicial;
	double* util = new double[this->numVertices];

	int* solucaoSucessora;

	while (k < this->iteracoes) {
		if (k % 5000 == 0) { // zerar penalidades
			for (int i = 0; i < this->numVertices - 1; i++)
				for (int j = 0; j < this->numVertices; j++)
					this->penalidade[i][j] = 0;
		}
		// atualizar matriz de sucessores
		switch (this->opcao) {
		case TSPsolverOpcao::OpcaoBuscaLocalRapidaArestasAleatorias:
			this->calcularMatrizSucessoresAleatoria();
			break;
		case TSPsolverOpcao::OpcaoBuscaLocalRapidaArestasMenoresPrimeiro:
			this->calcularMatrizSucessoresOrdenada(true);
			break;
		case TSPsolverOpcao::OpcaoBuscaLocalRapidaArestasMaioresPrimeiro:
			this->calcularMatrizSucessoresOrdenada(false);
			break;
		}

		std::printf("Iteração: %d\n", k);

		solucaoSucessora = this->fastLocalSearch(solucaoAtual, melhorSolucao, TSPsolverOpcao::OpcaoPrimeiroAprimorante,
				true);

		double maxUtil = 0;
		for (int i = 0; i < this->numVertices; i++) {
			util[i] = this->expressaoUtilidade(solucaoSucessora, i);
			if (util[i] > maxUtil)
				maxUtil = util[i];
		}
		for (int i = 0; i < this->numVertices; i++) {
			if (util[i] == maxUtil) {
				int atual = solucaoSucessora[i];
				int sucessor = solucaoSucessora[(i + 1) % this->numVertices];
				ajustaIndices(&atual, &sucessor);
				this->penalidade[atual][sucessor]++;

				// ativar sub-vizinhanças relacionadas aos vértices da aresta penalizada
				if (!this->featureSubNeighborhoodActivation[atual]) {
					this->featureSubNeighborhoodActivation[atual] = true;
					numSubNeighborhoodsActive++;
				}
				if (!this->featureSubNeighborhoodActivation[sucessor]) {
					this->featureSubNeighborhoodActivation[sucessor] = true;
					numSubNeighborhoodsActive++;
				}

				std::printf("--> aresta (%d, %d) penalizada: %d\n", atual + 1, sucessor + 1, this->penalidade[atual][sucessor]);
			}
		}

		k++;
	}

	delete[] solucaoAtual;
	delete[] solucaoSucessora;

	return melhorSolucao;
}

void efetuar2opt(int* rota, int* novaRota, int numVertices, int i, int j) {
	for (int i1 = 0; i1 < i; i1++) {
		novaRota[i1] = rota[i1];
	}
	for (int i2 = 0; i2 <= j - i; i2++) {
		novaRota[j - i2] = rota[i + i2];
	}
	for (int i3 = j + 1; i3 < numVertices; i3++) {
		novaRota[i3] = rota[i3];
	}
//	std::printf("  -> nova rota: ( ");
//	for (int k = 0; k < numVertices; k++) {
//		std::printf("%s%d ", (k > 0 ? "-> " : ""), novaRota[k] + 1);
//	}
//	std::printf(")\n");
}

// 2-opt
int* TSPsolver::localSearch(int* solucaoAtual, int* melhorSolucao, bool usarFuncaoCustoAumentada) {
	std::printf(" -> Efetuando Busca local...\n");
	int* otimoLocal = new int[this->numVertices];
	copiarVetor(solucaoAtual, otimoLocal, this->numVertices);
	int* melhorVizinho = new int[this->numVertices];
	copiarVetor(solucaoAtual, melhorVizinho, this->numVertices);

	double custoAumentadoOtimoLocal = this->funcaoCustoSolucaoAumentada(solucaoAtual);
	double custoOtimoLocal = this->funcaoCustoSolucao(solucaoAtual);
	double melhorCusto = this->funcaoCustoSolucao(melhorSolucao);
	double custoMelhorVizinho = custoOtimoLocal;
	double custoAumentadoMelhorVizinho = custoAumentadoOtimoLocal;

//	std::printf("  - melhor custo aumentado: %.1f\n    ...", custoAumentadoSolucaoAtual);

	int* novaSolucao = new int[this->numVertices];

	bool houveMelhora;
	int cont = 0;
	do {
		houveMelhora = false;
		for (int i = 0; i < this->numVertices; i++) {
			for (int j = i + 1; j < this->numVertices; j++) {
				if (i == 0 && j == this->numVertices - 1)
					continue; // 2-opt vai gerar solução equivalente à inicial...
				efetuar2opt(otimoLocal, novaSolucao, this->numVertices, i, j);
//				double custo = this->funcaoCustoSolucao(novaSolucao);
//				double custoAumentado = this->funcaoCustoSolucaoAumentada(novaSolucao);
				int antecessorI = (i > 0 ? i - 1 : this->numVertices - 1);
				int sucessorJ = (j + 1) % this->numVertices;
				int arestaRemovida1Origem = otimoLocal[antecessorI];
				int arestaRemovida1Destino = otimoLocal[i];
				int arestaRemovida2Origem = otimoLocal[j];
				int arestaRemovida2Destino = otimoLocal[sucessorJ];
				int arestaIncluida1Origem = arestaRemovida1Origem;
				int arestaIncluida1Destino = arestaRemovida2Origem;
				int arestaIncluida2Origem = arestaRemovida1Destino;
				int arestaIncluida2Destino = arestaRemovida2Destino;
				double custoArestaRemovida1 = this->obterPesoAresta(arestaRemovida1Origem, arestaRemovida1Destino);
				double custoArestaRemovida2 = this->obterPesoAresta(arestaRemovida2Origem, arestaRemovida2Destino);
				double custoArestaIncluida1 = this->obterPesoAresta(arestaIncluida1Origem, arestaIncluida1Destino);
				double custoArestaIncluida2 = this->obterPesoAresta(arestaIncluida2Origem, arestaIncluida2Destino);
				double custo = custoOtimoLocal - custoArestaRemovida1 - custoArestaRemovida2 + custoArestaIncluida1
						+ custoArestaIncluida2;
				double custoAumentado = custoAumentadoOtimoLocal - custoArestaRemovida1
						- this->lambda * this->obterPenalidadeAresta(arestaRemovida1Origem, arestaRemovida1Destino)
						- custoArestaRemovida2
						- this->lambda * this->obterPenalidadeAresta(arestaRemovida2Origem, arestaRemovida2Destino)
						+ custoArestaIncluida1
						+ this->lambda * this->obterPenalidadeAresta(arestaIncluida1Origem, arestaIncluida1Destino)
						+ custoArestaIncluida2
						+ this->lambda * this->obterPenalidadeAresta(arestaIncluida2Origem, arestaIncluida2Destino);

				if (custoAumentado < custoAumentadoMelhorVizinho) {
					copiarVetor(novaSolucao, melhorVizinho, this->numVertices);
					custoAumentadoMelhorVizinho = custoAumentado;
					custoMelhorVizinho = custo;
					houveMelhora = true;

//					std::printf("  *** melhoria no ótimo local, custo: %.1f\n", melhorCustoAumentado);
				}
				if (custo < melhorCusto) {
					copiarVetor(novaSolucao, melhorSolucao, this->numVertices);
					melhorCusto = custo;

//					std::printf("\n  *** melhoria no ótimo global, custo: %.1f\n    ...", melhorCusto);
				}
				cont++;
			}
		}
		if (houveMelhora) {
			copiarVetor(melhorVizinho, otimoLocal, this->numVertices);
			custoOtimoLocal = custoMelhorVizinho;
			custoAumentadoOtimoLocal = custoAumentadoMelhorVizinho;
		}
	} while (houveMelhora);

	std::printf("\n  - custo do ótimo local..........: %.1f\n", this->funcaoCustoSolucao(otimoLocal));
	std::printf("  - custo aumentado do ótimo local: %.1f\n", custoAumentadoOtimoLocal);
	std::printf("  - custo da melhor solução.......: %.1f\n", melhorCusto);
	std::printf(" -> Fim de busca local.\n");
	std::printf(" -> Soluções avaliadas: %d\n", cont);

	return otimoLocal;
}

double arredondar(double num, int casasDecimais) {
	int b = 1;
	for (int i = 0; i < casasDecimais; i++)
		b = b * 10;
	return floor(num * b + 0.5) / b;
}

int* TSPsolver::fastLocalSearch(int* solucaoAtual, int* melhorSolucao, int opcaoAprimorante,
		bool usarFuncaoCustoAumentada) {
	std::printf(" -> Efetuando Busca local...\n");
	int* otimoLocal = new int[this->numVertices];
	copiarVetor(solucaoAtual, otimoLocal, this->numVertices);

	double custoAumentadoOtimoLocal = this->funcaoCustoSolucaoAumentada(solucaoAtual);
	double custoOtimoLocal = this->funcaoCustoSolucao(solucaoAtual);
	double melhorCusto = this->funcaoCustoSolucao(melhorSolucao);

//	std::printf("  - melhor custo aumentado: %.1f\n    ...", custoAumentadoSolucaoAtual);

	int* novaSolucao = new int[this->numVertices];

	int* melhorVizinho = new int[this->numVertices];
	copiarVetor(solucaoAtual, melhorVizinho, this->numVertices);
	double custoAumentadoMelhorVizinho = custoAumentadoOtimoLocal;
	double custoMelhorVizinho = custoOtimoLocal;

	int cont = 0;

	int i = -1;
	while (numSubNeighborhoodsActive > 0) {
		i = (i + 1) % this->numVertices;
		int cidadeAtualPercurso = otimoLocal[i];

		if (this->featureSubNeighborhoodActivation[cidadeAtualPercurso]) {
			bool houveMelhora = false;
			for (int j = 0; j < this->numVertices; j++) {
				int indCidadeAtual = i;
				int indCidadeParaTrocar = this->sucessor[i][j];
				int cidadeParaTrocar = otimoLocal[indCidadeParaTrocar];
				if (indCidadeParaTrocar == indCidadeAtual)
					continue;
				if ((indCidadeAtual == 0 && indCidadeParaTrocar == this->numVertices - 1) || (indCidadeAtual == this->numVertices - 1 && indCidadeParaTrocar == 0))
					continue; // 2-opt vai gerar solução equivalente à inicial...
				ajustaIndices(&indCidadeAtual, &indCidadeParaTrocar);
				efetuar2opt(otimoLocal, novaSolucao, this->numVertices, indCidadeAtual, indCidadeParaTrocar);
				int antecessorCidadeAtual = (indCidadeAtual > 0 ? indCidadeAtual - 1 : this->numVertices - 1);
				int sucessorCidadeParaTrocar = (indCidadeParaTrocar + 1) % this->numVertices;
				int arestaRemovida1Origem = otimoLocal[antecessorCidadeAtual];
				int arestaRemovida1Destino = otimoLocal[indCidadeAtual];
				int arestaRemovida2Origem = otimoLocal[indCidadeParaTrocar];
				int arestaRemovida2Destino = otimoLocal[sucessorCidadeParaTrocar];
				int arestaIncluida1Origem = arestaRemovida1Origem;
				int arestaIncluida1Destino = arestaRemovida2Origem;
				int arestaIncluida2Origem = arestaRemovida1Destino;
				int arestaIncluida2Destino = arestaRemovida2Destino;
				double custoArestaRemovida1 = this->obterPesoAresta(arestaRemovida1Origem, arestaRemovida1Destino);
				double custoArestaRemovida2 = this->obterPesoAresta(arestaRemovida2Origem, arestaRemovida2Destino);
				double custoArestaIncluida1 = this->obterPesoAresta(arestaIncluida1Origem, arestaIncluida1Destino);
				double custoArestaIncluida2 = this->obterPesoAresta(arestaIncluida2Origem, arestaIncluida2Destino);
				double custo = custoOtimoLocal - custoArestaRemovida1 - custoArestaRemovida2 + custoArestaIncluida1
						+ custoArestaIncluida2;
				double custoAumentado = custoAumentadoOtimoLocal - custoArestaRemovida1
						- this->lambda * this->obterPenalidadeAresta(arestaRemovida1Origem, arestaRemovida1Destino)
						- custoArestaRemovida2
						- this->lambda * this->obterPenalidadeAresta(arestaRemovida2Origem, arestaRemovida2Destino)
						+ custoArestaIncluida1
						+ this->lambda * this->obterPenalidadeAresta(arestaIncluida1Origem, arestaIncluida1Destino)
						+ custoArestaIncluida2
						+ this->lambda * this->obterPenalidadeAresta(arestaIncluida2Origem, arestaIncluida2Destino);

				if ((opcaoAprimorante == TSPsolverOpcao::OpcaoPrimeiroAprimorante &&
						arredondar(custoAumentado, 4) < arredondar(custoAumentadoOtimoLocal, 4)) ||
						(opcaoAprimorante == TSPsolverOpcao::OpcaoMelhorAprimorante &&
						arredondar(custoAumentado, 4) < arredondar(custoAumentadoMelhorVizinho, 4))) {
					if (opcaoAprimorante == TSPsolverOpcao::OpcaoPrimeiroAprimorante) {
						copiarVetor(novaSolucao, otimoLocal, this->numVertices);
						custoAumentadoOtimoLocal = custoAumentado;
						custoOtimoLocal = custo;
					} else if (opcaoAprimorante == TSPsolverOpcao::OpcaoMelhorAprimorante) {
						copiarVetor(novaSolucao, melhorVizinho, this->numVertices);
						custoAumentadoMelhorVizinho = custoAumentado;
						custoMelhorVizinho = custo;
					}

					houveMelhora = true;

					// ativar sub-vizinhança
					if (!this->featureSubNeighborhoodActivation[cidadeParaTrocar]) {
						this->featureSubNeighborhoodActivation[cidadeParaTrocar] = true;
						numSubNeighborhoodsActive++;
					}
				}
				if (custo < melhorCusto) {
					copiarVetor(novaSolucao, melhorSolucao, this->numVertices);
					melhorCusto = custo;
				}
				cont++;

				if (houveMelhora && opcaoAprimorante == TSPsolverOpcao::OpcaoPrimeiroAprimorante) {
					// atualizar ótimo local
					copiarVetor(novaSolucao, otimoLocal, this->numVertices);
					custoAumentadoOtimoLocal = custoAumentado;
					custoOtimoLocal = custo;

					// backtrack one city
					i--;

					// sair da busca por esta vizinhança para retomar a busca na próxima iteração do while
					break;
				}

			}
			if (!houveMelhora) {
				// desativar sub-vizinhança
				this->featureSubNeighborhoodActivation[cidadeAtualPercurso] = false;
				numSubNeighborhoodsActive--;
			} else if (opcaoAprimorante == TSPsolverOpcao::OpcaoMelhorAprimorante) {
				// atualizar ótimo local
				copiarVetor(melhorVizinho, otimoLocal, this->numVertices);
				custoAumentadoOtimoLocal = custoAumentadoMelhorVizinho;
				custoOtimoLocal = custoMelhorVizinho;
			}
		}
	}

	std::printf("\n  - custo do ótimo local..........: %.1f\n", this->funcaoCustoSolucao(otimoLocal));
	std::printf("  - custo aumentado do ótimo local: %.1f\n", custoAumentadoOtimoLocal);
	std::printf("  - custo da melhor solução.......: %.1f\n", melhorCusto);
	std::printf(" -> Fim de busca local.\n");
	std::printf(" -> Soluções avaliadas: %d\n", cont);

	return otimoLocal;
}

double TSPsolver::obterPesoAresta(int origem, int destino) {
	return this->pesos[std::min(origem, destino)][std::max(origem, destino)];
}

int TSPsolver::obterPenalidadeAresta(int origem, int destino) {
	return this->penalidade[std::min(origem, destino)][std::max(origem, destino)];
}

int* TSPsolver::gerarSolucaoArbitraria() {
	srand(time(NULL));
	int* solucao = new int[this->numVertices];
	for (int i = 0; i < this->numVertices; i++)
		solucao[i] = -1;
	for (int i = 0; i < this->numVertices; i++) {
		int r = rand() % this->numVertices;
		if (solucao[r] == -1) // se posição r estiver livre, ocupe com i
			solucao[r] = i;
		else { // senão, busque a próxima posição livre a partir de r e ocupe com i
			for (int j = 1; j < this->numVertices; j++) {
				if (solucao[(r + j) % this->numVertices] == -1) {
					solucao[(r + j) % this->numVertices] = i;
					break;
				}
			}
		}
	}

	return solucao;
}
