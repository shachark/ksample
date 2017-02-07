/** ============================================================================

    Description: Dynamic Slicing for Hypothesis Testing
         Author: Chao Ye, Bioinformatics Division, TNLIST, Tsinghua University
                 Bo Jiang, Department of Statistics, Harvard University
          Email: yec09@mails.tsinghua.edu.cn, bjiang@fas.harvard.edu
    Last Update: Feb 14, 2014

============================================================================ **/

#include "Dynslicing.h"

using namespace std;

/**
  One-sample hypothesis testing
  Dynamic slicing based on equipartition (one over n resolution)
  y is the quantile according to F_0, y \in (0, 1)
**/
double dynslicing_1eqp(double *y, int len, double lambda)
{
	double lpd = -lambda * log(len);
	double unit = 1.0 / len;
	int *ctab = new int[len+1];  //  cumulated unit counts
	ctab[0] = 0;
	int l = 1;
	double flagf = unit;
	for(int i = 0; i < len; ++i){
		while(y[i] > flagf){
			ctab[l++] = i;
			flagf += unit;
		}
	}
	for(; l < len+1; ++l){
		ctab[l] = len;
	}
	double *score = new double[len+1];
	int *idx = new int[len+1];
	for(int k = 0; k < len+1; ++k){
		score[k] = 0;
		idx[k] = -1;
	}
	int cutpos;
	double dq, dc, tpcut, cutsc;
	for(int i = 1; i < len+1; ++i){
		cutsc = lpd + score[0];
		dq = i;
		dc = ctab[i] - ctab[0];
		if(dc > EPSILON){
			cutsc +=  dc * log(dc / dq);
		}
		cutpos = 0;
		for(int j = 1; j < i; ++j){
			tpcut = lpd + score[j];
			dq = i - j;
			dc = ctab[i] - ctab[j];
			if(dc > EPSILON){
				tpcut +=  dc * log(dc / dq);
			}
			if(cutsc < tpcut){
				cutsc = tpcut;
				cutpos = j;
			}
		}
		score[i] = cutsc;
		idx[i] = cutpos;
	}
	double mlik = score[len] - lpd;
	dq = len;
	dc = ctab[len];
	mlik -= dc * log(dc / dq);
	delete[] ctab;
	delete[] score;
	delete[] idx;
	return mlik;  //  maximized log-likelihood
}

/**
  Dynamic slicing for one-sample hypothesis testing
  y is the quantile according to F_0, y \in (0, 1)
**/
double dynslicing_one(double *y, int len, double lambda, double alpha)
{
	double lpd = -lambda * log(len);
	double unit = 1.0 / len;
	double *lscore = new double[len+1];
	double *rscore = new double[len+1];
	int *idx = new int[2*len+2];
	for(int k = 0; k < len+1; ++k){
		lscore[k] = 0;
		rscore[k] = 0;
		idx[k] = -1;
	}
	int cutpos;
	double dc, dq, lval, rval;
	double tpcut, cutsc;
	//  i = 0, log likelihood = 0 for the left side
	lval = alpha * log(y[0]);
	rval = log(unit) + (alpha - 1) * log(y[0]);	//	log(unit / y[0]) + alpha * log(y[0])
	lscore[0] = lval + lpd;
	rscore[0] = rval + lpd;
	idx[1] = 0;
	idx[2] = 0;
	//  i = 0 end
	for(int i = 1; i < len; ++i){
		//  on left side of i
		lval = i * log(i * unit / y[i]);
		cutsc = lval + lpd + alpha * log(y[i]);
		cutpos = 0;
		for(int j = 0; j < i-1; ++j){
			dc = i - j;
			dq = y[i] - y[j];
			lval = lscore[j] + dc * log(dc * unit / dq);
			rval = rscore[j] + (dc - 1) * log((dc - 1) * unit / dq);
			tpcut = max(lval, rval) + lpd + alpha * log(dq);
			if(cutsc < tpcut){
				cutsc = tpcut;
				cutpos = j;
			}
		}
		//  j = i - 1
		dc = 1;
		dq = y[i] - y[i-1];
		lval = lscore[i-1] + log(unit / dq);
		rval = rscore[i-1];
		tpcut = max(lval, rval) + lpd + alpha * log(dq);
		if(cutsc < tpcut){
			cutsc = tpcut;
			cutpos = i-1;
		}
		//  j = i - 1 end
		lscore[i] = cutsc;
		idx[2*i+1] = cutpos;
		//  on right side of i
		rval = (i + 1) * log((i + 1) * unit / y[i]);
		cutsc = rval + lpd + alpha * log(y[i]);
		cutpos = 0;
		for(int j = 0; j < i; ++j){
			dc = i - j;
			dq = y[i] - y[j];
			lval = lscore[j] + (dc + 1) * log((dc + 1) * unit / dq);
			rval = rscore[j] + dc * log(dc * unit / dq);
			tpcut = max(lval, rval) + lpd + alpha * log(dq);
			if(cutsc < tpcut){
				cutsc = tpcut;
				cutpos = j;
			}
		}
		rscore[i] = cutsc;
		idx[2*i+2] = cutpos;
	}
	//  i = len, i.e. for y_{n+1} = 1.0
	//  on left side of len, no right side
	cutsc = lpd;  //  n*log(n/n /1) + lpd + alpha * log(1)
	cutpos = 0;
	for(int j = 0; j < len-1; ++j){
		dc = len - j;
		dq = 1.0 - y[j];
		lval = lscore[j] + dc * log(dc * unit / dq);
		rval = rscore[j] + (dc - 1) * log((dc - 1) * unit / dq);
		tpcut = max(lval, rval) + lpd + alpha * log(dq);
		if(cutsc < tpcut){
			cutsc = tpcut;
			cutpos = j;
		}
	}
	//  j = len - 1
	dc = 1;
	dq = 1.0 - y[len-1];
	lval = lscore[len-1] + log(unit / dq);
	rval = rscore[len-1];
	tpcut = max(lval, rval) + lpd + alpha * log(dq);
	if(cutsc < tpcut){
		cutsc = tpcut;
		cutpos = len-1;
	}
	//  j = len - 1 end
	lscore[len] = cutsc;
	idx[2*len+1] = cutpos;
	//  i = len end
	double mlik = lscore[len] - lpd;
	delete[] lscore;
	delete[] rscore;
	delete[] idx;

	return mlik;
}

/**
  K-sample hypothesis testing
  Dynamic slicing based on equipartition (one over sqrt{n} resolution)
  X: discrete variable; Y: continuous variable.
  X is sorted according to value of Y in advanced.
  Never slice within a clump when maximizing the regularized likelihood-ratio.
**/
double dynslicing_keqp(int *x, int len, int dim, double lambda)
{
	double lpd = -lambda * log(len);  //  penalty (log of prior odds) for each additional slice
	int baselen = sqrt(len);
	int *divsch = new int[baselen + 3];  //  at most "baselen+1" groups with each group has size "baselen" at least
	divsch[0] = 0;
	int flagl = baselen;
	int ngrp = 0;
	while(flagl < len){
		if(x[flagl] - x[flagl-1] != 0){
			divsch[++ngrp] = flagl;
			flagl += baselen;
		}else{
			flagl++;
		}
	}
	divsch[++ngrp] = len;
	int **ctab = new int*[ngrp+1];
	for(int k = 0; k < ngrp+1; ++k){
		ctab[k] = new int[dim];
		for(int j = 0; j < dim; ++j){
			ctab[k][j] = 0;
		}
	}
	for(int k = 1; k < ngrp+1; ++k){
		for(int j = 0; j < dim; ++j){
			ctab[k][j] += ctab[k-1][j];
		}
		for(int j = divsch[k-1]; j < divsch[k]; ++j){
			ctab[k][x[j]] += 1;
		}
	}
	double *score = new double[ngrp+1];
	int *idx = new int[ngrp+1];
	for(long k = 0; k < ngrp+1; ++k){
		score[k] = 0;
		idx[k] = -1;
	}
	double *counts = new double[dim];
	int tc, cutpos;
	double tpcut, cutsc;
	for(int i = 1; i < ngrp+1; ++i){
		//  j = 0
		cutsc = lpd + score[0];
		tc = 0;
		for(int u = 0; u < dim; ++u){
			tc += ctab[i][u];
		}
		for(int u = 0; u < dim; ++u){
			counts[u] = (double)(ctab[i][u]);
			if(counts[u] > EPSILON){
				cutsc += counts[u] * log(counts[u] / tc);
			}
		}
		cutpos = 0;
		//  j = 0 end
		for(int j = 1; j < i; ++j){
			tpcut = lpd + score[j];
			tc = 0;
			for(int u = 0; u < dim; ++u){
				tc += ctab[i][u] - ctab[j][u];
			}
			for(int u = 0; u < dim; ++u){
				counts[u] = (double)(ctab[i][u] - ctab[j][u]);
				if(counts[u] > EPSILON){
					tpcut += counts[u] * log(counts[u] / tc);
				}
			}
			if(cutsc < tpcut){
				cutsc = tpcut;
				cutpos = j;
			}
		}
		score[i] = cutsc;
		idx[i] = cutpos;
	}
	delete[] counts;
	double mlik = score[ngrp] - lpd;
	//  substract null log-likelihood (assume one slice, i.e., no cut)
	for(int u = 0; u < dim; ++u){
		mlik -= ctab[ngrp][u] * log((double)(ctab[ngrp][u]) / len);
	}
	delete[] score;
	delete[] idx;
	for(int k = 0; k < ngrp; ++k){
		delete ctab[k];
	}
	delete[] ctab;
	delete[] divsch;
	return mlik;  //  maximized log-likelihood
}

double dynslicing_keqp(int *x, int len, int dim, double lambda, vector<int> &slices)
{
	double lpd = -lambda * log(len);  //  penalty (log of prior odds) for each additional slice
	int baselen = sqrt(len);
	int *divsch = new int[baselen + 3];  //  at most "baselen+1" groups with each group has size "baselen" at least
	divsch[0] = 0;
	int flagl = baselen;
	int ngrp = 0;
	while(flagl < len){
		if(x[flagl] - x[flagl-1] != 0){
			divsch[++ngrp] = flagl;
			flagl += baselen;
		}else{
			flagl++;
		}
	}
	divsch[++ngrp] = len;
	int **ctab = new int*[ngrp+1];
	for(int k = 0; k < ngrp+1; ++k){
		ctab[k] = new int[dim];
		for(int j = 0; j < dim; ++j){
			ctab[k][j] = 0;
		}
	}
	for(int k = 1; k < ngrp+1; ++k){
		for(int j = 0; j < dim; ++j){
			ctab[k][j] += ctab[k-1][j];
		}
		for(int j = divsch[k-1]; j < divsch[k]; ++j){
			ctab[k][x[j]] += 1;
		}
	}
	double *score = new double[ngrp+1];
	int *idx = new int[ngrp+1];
	for(long k = 0; k < ngrp+1; ++k){
		score[k] = 0;
		idx[k] = -1;
	}
	double *counts = new double[dim];
	int tc, cutpos;
	double tpcut, cutsc;
	for(int i = 1; i < ngrp+1; ++i){
		//  j = 0
		cutsc = lpd + score[0];
		tc = 0;
		for(int u = 0; u < dim; ++u){
			tc += ctab[i][u];
		}
		for(int u = 0; u < dim; ++u){
			counts[u] = (double)(ctab[i][u]);
			if(counts[u] > EPSILON){
				cutsc += counts[u] * log(counts[u] / tc);
			}
		}
		cutpos = 0;
		//  j = 0 end
		for(int j = 1; j < i; ++j){
			tpcut = lpd + score[j];
			tc = 0;
			for(int u = 0; u < dim; ++u){
				tc += ctab[i][u] - ctab[j][u];
			}
			for(int u = 0; u < dim; ++u){
				counts[u] = (double)(ctab[i][u] - ctab[j][u]);
				if(counts[u] > EPSILON){
					tpcut += counts[u] * log(counts[u] / tc);
				}
			}
			if(cutsc < tpcut){
				cutsc = tpcut;
				cutpos = j;
			}
		}
		score[i] = cutsc;
		idx[i] = cutpos;
	}
	delete[] counts;
	
	int slicenum = 0;
	int flag = ngrp;
	while(flag > 0){
		flag = idx[flag];
		slicenum++;
	}
	vector<int>().swap(slices);
	slices.resize((dim + 1) * slicenum);
	flag = ngrp;
	vector<int> spos(slicenum+1);
	for(int i = slicenum; i > -1; --i){
		spos[i] = flag;
		flag = idx[flag];
	}
	spos[0] = 0;
	int subcounts = 0;
	for(int i = 0; i < slicenum; ++i){
		subcounts = 0;
		for(int j = 0; j < dim; ++j){
			slices[i*(dim+1) + j] = ctab[spos[i+1]][j] - ctab[spos[i]][j];
			subcounts += slices[i*(dim+1) + j];
		}
		slices[i*(dim+1) + dim] = subcounts;
	}

	double mlik = score[ngrp] - lpd;
	//  substract null log-likelihood (assume one slice, i.e., no cut)
	for(int u = 0; u < dim; ++u){
		mlik -= ctab[ngrp][u] * log((double)(ctab[ngrp][u]) / len);
	}
	delete[] score;
	delete[] idx;
	for(int k = 0; k < ngrp; ++k){
		delete ctab[k];
	}
	delete[] ctab;
	delete[] divsch;
	return mlik;  //  maximized log-likelihood
}

/**
  Dynamic slicing for K-sample hypothesis testing
  X: discrete variable; Y: continuous variable.
  X is sorted according to value of Y in advanced.
  Never slice within a clump when maximizing the regularized likelihood-ratio.
**/
double dynslicing_k(int *x, int len, int dim, double lambda)
{
	double lpd = -lambda * log(len);  //  penalty (log of prior odds) for each additional slice
	int **ctab = new int*[len+1];
	for(int k = 0; k < len+1; ++k){
		ctab[k] = new int[dim];
		for(int j = 0; j < dim; ++j){
			ctab[k][j] = 0;
		}
	}
	int flagl = 1;
	int clpcount = 1;  //  clump count
	int clumpnum = 1;  //  clump number
	while(flagl < len){
		if(x[flagl] - x[flagl-1] != 0){
			ctab[clumpnum][x[flagl-1]] = clpcount;
			clpcount = 1;
			clumpnum++;
		}else{
			clpcount++;
		}
		flagl++;
	}
	ctab[clumpnum][x[len-1]] = clpcount;

	for(int k = 1; k < clumpnum+1; ++k){
		for(int j = 0; j < dim; ++j){
			ctab[k][j] += ctab[k-1][j];
		}
	}
	//  dynamic programming
	double *score = new double[clumpnum+1];
	int *idx = new int[clumpnum+1];
	for(int k = 0; k < clumpnum+1; ++k){
		score[k] = 0;
		idx[k] = -1;
	}
	double *counts = new double[dim];
	int tc, cutpos;
	double tpcut, cutsc;
	for(int i = 1; i < clumpnum+1; ++i){
		//  j = 0
		cutsc = lpd + score[0];
		tc = 0;
		for(int u = 0; u < dim; ++u){
			tc += ctab[i][u];
		}
		for(int u = 0; u < dim; ++u){
			counts[u] = (double)(ctab[i][u]);
			if(counts[u] > EPSILON){
				cutsc += counts[u] * log(counts[u] / tc);
			}
		}
		cutpos = 0;
		//  j = 0 end
		for(int j = 1; j < i; ++j){
			tpcut = lpd + score[j];
			tc = 0;
			for(int u = 0; u < dim; ++u){
				tc += ctab[i][u] - ctab[j][u];
			}
			for(int u = 0; u < dim; ++u){
				counts[u] = (double)(ctab[i][u] - ctab[j][u]);
				if(counts[u] > EPSILON){
					tpcut += counts[u] * log(counts[u] / tc);
				}
			}
			if(cutsc < tpcut){
				cutsc = tpcut;
				cutpos = j;
			}
		}
		score[i] = cutsc;
		idx[i] = cutpos;
	}
	double mlik = score[clumpnum] - lpd;
	//  substract null log-likelihood (assume one slice, i.e., no cut)
	for(int u = 0; u < dim; ++u){
		if(ctab[clumpnum][u] > EPSILON){
			mlik -= ctab[clumpnum][u] * log((double)(ctab[clumpnum][u]) / len);
		}
	}
	delete[] score;
	delete[] idx;
	delete[] counts;
	for(int k = 0; k < len; ++k){
		delete ctab[k];
	}
	delete[] ctab;
	return mlik;  //  maximized log-likelihood
}

double dynslicing_k(int *x, int len, int dim, double lambda, vector<int> &slices)
{
	double lpd = -lambda * log(len);  //  penalty (log of prior odds) for each additional slice
	int **ctab = new int*[len+1];
	for(int k = 0; k < len+1; ++k){
		ctab[k] = new int[dim];
		for(int j = 0; j < dim; ++j){
			ctab[k][j] = 0;
		}
	}
	int flagl = 1;
	int clpcount = 1;  //  clump count
	int clumpnum = 1;  //  clump number
	while(flagl < len){
		if(x[flagl] - x[flagl-1] != 0){
			ctab[clumpnum][x[flagl-1]] = clpcount;
			clpcount = 1;
			clumpnum++;
		}else{
			clpcount++;
		}
		flagl++;
	}
	ctab[clumpnum][x[len-1]] = clpcount;

	for(int k = 1; k < clumpnum+1; ++k){
		for(int j = 0; j < dim; ++j){
			ctab[k][j] += ctab[k-1][j];
		}
	}
	//  dynamic programming
	double *score = new double[clumpnum+1];
	int *idx = new int[clumpnum+1];
	for(int k = 0; k < clumpnum+1; ++k){
		score[k] = 0;
		idx[k] = -1;
	}
	double *counts = new double[dim];
	int tc, cutpos;
	double tpcut, cutsc;
	for(int i = 1; i < clumpnum+1; ++i){
		//  j = 0
		cutsc = lpd + score[0];
		tc = 0;
		for(int u = 0; u < dim; ++u){
			tc += ctab[i][u];
		}
		for(int u = 0; u < dim; ++u){
			counts[u] = (double)(ctab[i][u]);
			if(counts[u] > EPSILON){
				cutsc += counts[u] * log(counts[u] / tc);
			}
		}
		cutpos = 0;
		//  j = 0 end
		for(int j = 1; j < i; ++j){
			tpcut = lpd + score[j];
			tc = 0;
			for(int u = 0; u < dim; ++u){
				tc += ctab[i][u] - ctab[j][u];
			}
			for(int u = 0; u < dim; ++u){
				counts[u] = (double)(ctab[i][u] - ctab[j][u]);
				if(counts[u] > EPSILON){
					tpcut += counts[u] * log(counts[u] / tc);
				}
			}
			if(cutsc < tpcut){
				cutsc = tpcut;
				cutpos = j;
			}
		}
		score[i] = cutsc;
		idx[i] = cutpos;
	}
	
	int slicenum = 0;
	int flag = clumpnum;
	while(flag > 0){
		flag = idx[flag];
		slicenum++;
	}
	vector<int>().swap(slices);
	slices.resize((dim + 1) * slicenum);
	flag = clumpnum;
	vector<int> spos(slicenum+1);
	for(int i = slicenum; i > -1; --i){
		spos[i] = flag;
		flag = idx[flag];
	}
	spos[0] = 0;
	int subcounts = 0;
	for(int i = 0; i < slicenum; ++i){
		subcounts = 0;
		for(int j = 0; j < dim; ++j){
			slices[i*(dim+1) + j] = ctab[spos[i+1]][j] - ctab[spos[i]][j];
			subcounts += slices[i*(dim+1) + j];
		}
		slices[i*(dim+1) + dim] = subcounts;
	}
	
	double mlik = score[clumpnum] - lpd;
	//  substract null log-likelihood (assume one slice, i.e., no cut)
	for(int u = 0; u < dim; ++u){
		if(ctab[clumpnum][u] > EPSILON){
			mlik -= ctab[clumpnum][u] * log((double)(ctab[clumpnum][u]) / len);
		}
	}
	delete[] score;
	delete[] idx;
	delete[] counts;
	for(int k = 0; k < len; ++k){
		delete ctab[k];
	}
	delete[] ctab;
	return mlik;  //  maximized log-likelihood
}
