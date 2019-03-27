//Copyright(c) 2015 M. Mancini
//
//Permission is hereby granted, free of charge, to any person obtaining a copy
//of this software and associated documentation files(the "Software"), to deal
//in the Software without restriction, including without limitation the rights
//to use, copy, modify, merge, publish, distribute, sublicense, and / or sell
//copies of the Software, and to permit persons to whom the Software is
//furnished to do so, subject to the following conditions :
//
//The above copyright notice and this permission notice shall be included in all
//copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE
//AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//SOFTWARE.

#pragma once

class StableFit
{
public:
	StableFit();
	~StableFit();
	
	std::vector<double> Execute(std::vector<double> data);
	std::vector<double> intAlphaBeta(std::vector<double> data);
	std::vector<double> intGamDel(std::vector<double> X, double alpha, double beta);
	
	double interp2(double *a, double *b, double *v, double ta, double tb, int rows, int cols);
	double interp3(double *a, double *b, double *c,  double *v, double ta, double tb, double tc, int rows, int cols, int depth);
	std::vector<double> prctile(std::vector<double> data, std::vector<double> prcvec);
	std::vector<double> regress(std::vector<double> y, std::vector<double> x, int nvars, int nsamples);
	void meshgrid(std::vector<double> xgv, std::vector<double> ygv, std::vector<double> zgv, double ***&X, double ***&Y, double ***&Z);
	std::vector<double> intGuess(double alpha, double beta, std::vector<double> u);
	std::vector<double> stblinv(std::vector<double> u, double alpha, double beta, double gam, double delta);
    
};

