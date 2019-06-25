#include "matplotlibcpp.h"
#include <iostream>
#include <stdlib.h>
#include <vector>
#include "string.h"
#include "windowing.cpp"
#include "Wavreader.cpp"
#include <fftw3.h>
#include "resample.h"
#include <Eigen/Core>

using namespace std;

int main(){

	//init variables
	bool avg = 1;

	//1.1 Read the audio file and sampling rate
	int sampleRate = readWav();
	int sampleRateTarget = 16000;
	cout << "Sample rate is: " << sampleRate << endl;

	vector<double> dataRaw = {};
	ifstream inFile;
	inFile.open("data/raw.dat");
	double tp;
	while(inFile >> tp){
		dataRaw.push_back(tp);
	}											//All data moved from list.dat to dataV, kinda unnecessary to put it in textfile, could fix
	inFile.close();

	//1.2 Make sure the sampling rate is 16kHz, resample if not
	vector<double> dataV = {};
	vector<double> dataResampled = {};
	if(sampleRate != sampleRateTarget){
		resample(sampleRateTarget,sampleRate,dataRaw,dataV);
	}
	else {
		for(auto i = dataRaw.begin(); i != dataRaw.end(); i++) dataV.push_back(*i);
	}

	//1.3 Split the data sequence into windows, implement windowing.cpp
	int frameLength = sampleRate * 0.03; 	//flattop : 0.03
	int hopSize 	= sampleRate * 0.02; 	//flattop : 0.02 (2/3 overlap)
	vector<vector<double>> frameMatrix = window("data/raw.dat", frameLength, hopSize, "flattop");

	//Plots
	namespace plt = matplotlibcpp;
	using namespace plt;

    figure(1);
	subplot(3,1,1);
	vector<double> t_axis1 = {};  //divide 0 - length of dataV by 16000
	int size = dataV.size();
	for(int i = 0; i < size; i++){
		t_axis1.push_back(i/sampleRateTarget);
	}
	plot(t_axis1, dataV);
 	xlim(0,2);


	subplot(3,1,2);
	vector<double> t_axis2 = {}; //divide 0 - frameLength with 16
	for(float i = 0; i < frameLength; i++){
		t_axis2.push_back(i/(sampleRate/1000));
	}
	vector<double> temp = {};
	unsigned int count1 = 0;
	for(int i = 0; i < frameLength; i++){
		temp.push_back(frameMatrix[22][i]); //frameMatrix[:,23]?
		count1++;
		if(count1 == t_axis2.size()) break;
	}
	plot(t_axis2, temp); 
	ylim(-1500, 2000);

	//1.4 FFTW
	double matrixSize = frameMatrix.size()*frameLength;					//total size of frameMatrix (n x m)
	double* in = (double*) fftw_malloc(sizeof(double) * matrixSize);	//allocate fft-input array
	fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * matrixSize/2 + 1);	//allocate fft-output array
    fftw_plan p;														//define fft plan
	int N = 1024;
	ofstream offtfile;
	offtfile.open("data/fft.dat"); 										//open textfile for data to matlab for testing
	vector<vector<double>> outV = {};									//vector for fft-output data
	vector<double> outVCol = {};
    p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_MEASURE);					//init fft-plan
	for(unsigned int i = 0; i < frameMatrix.size(); i++){				//for each column in frameMatrix, set it as fft-input
		for(int j = 0; j < frameLength; j++){							//
			in[j] = frameMatrix[i][j];									//
		}																//
		fftw_execute(p);												//execute fft for each column
		for(int i = N/2+1; i > 0; i--){									//
			offtfile << 20*log10(abs(out[i][0] + out[i][1])) << endl;	//send output to textfile
			outVCol.push_back(20*log10(abs(out[i][0] + out[i][1])));	//push output to outV - datavector
		}
		outV.push_back(outVCol);
		outVCol.clear();
	}   
	offtfile.close();		//close outputfile
	fftw_destroy_plan(p);	//destruct fft_plan
    fftw_free(in); 			//free allocated memory
	fftw_free(out);			//

	unsigned int nEnerBands = 32;
	vector<vector<double>> enerOut = {};
	vector<double> enerFrame = {};
	//apply ener_bands to outV
	int nfft = N/2 + 1;			//513
	int nframes = outV.size(); 	//216
	int bandsPerEner = floor(nfft/nEnerBands);
	unsigned int count = 0;

	for(int i = nfft-1; i > 0; i--){
		if(count % bandsPerEner == 0){
			int mean = 0;
			for(int j = 0; j < nframes; j++){
				int sum = 0;
				for (int k = 0; k < bandsPerEner; k++){
					if(i-k >= 0) sum += outV[j][i-k];
				} 
				sum = 10*log10(pow(abs(sum),2));
				mean += sum;										
				enerFrame.push_back(sum);													
			}	
			mean /= enerFrame.size();
			for(unsigned int k = 0; k < enerFrame.size(); k++){			
				enerFrame[k] -= mean;
			}	
			enerOut.push_back(enerFrame);
			enerFrame.clear();	
		}	
		count++;	
	}
	ofstream oEnerFile;
	oEnerFile.open("data/energy.dat");
	for(int i = 0; i < nframes; i++){
		for(int j = 0; j < 32; j++){
		oEnerFile << enerOut[j][i] << endl;	
		}
	}
	oEnerFile.close();

	//int avgLen;
	vector<vector<double>> avgEner = {};
	vector<double> avgEnerFrame = {};
	if(avg){
		//avgLen = ceil(enerOut[0].size()/5);
		for(unsigned int i = 0; i < enerOut[0].size(); i += 5){
			unsigned int idxEnd;
			if(i+5 < enerOut[0].size()) idxEnd = i+5;
			else idxEnd = enerOut[0].size();	
			for(unsigned int j = 0; j < enerOut.size(); j++){		
				int mean = 0;
				int div = 0;
				for(unsigned int k = i; k < idxEnd; k++){			
					mean += enerOut[j][k];
					div++;
				}
				if(div != 0) mean /= div;
				avgEnerFrame.push_back(mean); 
			}
			avgEner.push_back(avgEnerFrame);
			avgEnerFrame.clear();
		}
	}
	ifstream iTMatFile, iCtxOptFile, iCtxShapeFile;
	iTMatFile.open("data/T_mat.txt");
	iCtxOptFile.open("data/ctx_opt.txt");	
	iCtxShapeFile.open("data/ctx_shape.txt");	
		
	vector<double> ctxOptData, ctxShapeData = {};
	vector<Eigen::MatrixXd> tRealMatData;
	double tMatData[9][9][30] = {};
	double a, b, c;
	
	for(int i = 0; i < 9; i++){
		for(int j = 0; j < 30; j++){
			for(int k = 0; k < 9; k++){
				if(iTMatFile >> a){
					tMatData[k][i][j] = a;
				}	
			}
		}
	}
	
	Eigen::MatrixXd m(9,9);
	for(int i = 0; i < 30; i++){
		for(int j = 0; j < 9; j++){
			for(int k = 0; k < 9; k++){
				m(j, k)= tMatData[k][j][i];
			}
		}
	tRealMatData.push_back(m.transpose()); //tRealMatData  = T_mat = transformation matrix in correct order
	}
	
	while(iCtxOptFile >> b){
		ctxOptData.push_back(b);
	}
	while(iCtxShapeFile >> c){
		ctxShapeData.push_back(c);
	}
	iTMatFile.close();
	iCtxOptFile.close();
	iCtxShapeFile.close();

	//limits, first and last bands and first and last frames ignored
	int yMin = ctxShapeData[0];
	int yMax = ctxShapeData[0];
	int xMin = ctxShapeData[1];
	int xMax = ctxShapeData[1];
	for(unsigned int i = 0; i < ctxShapeData.size(); i += 2){
		if(i % 2 == 0){
			if(ctxShapeData[i] < yMin) yMin = ctxShapeData[i];
			else if(ctxShapeData[i] > yMax) yMax = ctxShapeData[i];
		}		
		else {
			if(ctxShapeData[i] < xMin) xMin = ctxShapeData[i];
			else if(ctxShapeData[i] > xMax) xMax = ctxShapeData[i];
		}	
	}
	yMin = -yMin;
	xMin = -xMin;

	int nEnerBandsUsed = nEnerBands - yMin - yMax;
	int nFramesUsed = avgEner.size()/*enerOut[0].size()*/ - xMin - xMax;

	double ctxBands[nFramesUsed][ctxShapeData.size()/2][nEnerBandsUsed] = {};
	double multAux [nFramesUsed][ctxShapeData.size()/2][nEnerBandsUsed] = {};
	double decorMat[nFramesUsed][nEnerBandsUsed] = {};

	Eigen::VectorXd v(9); 
	Eigen::MatrixXd x;

	for(int i = 0; i < nEnerBandsUsed; i++){								//30
		for(int j = 0; j < nFramesUsed-2; j++){ 							//25
			int mean = 0;	
			for(unsigned int k = 0; k < ctxShapeData.size()/2; k++){		//9
				int ctxIdx = ctxOptData[k*nEnerBandsUsed + j] - 1;
				int idx = j + xMin + ctxShapeData[ctxIdx*2 + 1];
				int idy = i + yMin + ctxShapeData[ctxIdx*2];
				ctxBands[j][k][i] = avgEner[idx][idy];//enerOut[idy][idx];
				v(k) = ctxBands[j][k][i];
				mean += v(k);
			}
			mean /= 9;
			for(int k = 0; k < 9; k++){
				v(k) -= mean;
			}
			v.transpose();
			x = tRealMatData[i]*v;
			x.transpose();
			for(unsigned int k = 0; k < ctxShapeData.size()/2; k++){
				multAux[j][k][i] = x(k);
			}
			decorMat[j][i] = multAux[j][0][i];
		}
	}
	ofstream outputDecorFile;
	outputDecorFile.open("data/decor.dat");
	cout << x.rows() << "x" << x.cols() << endl;
	for(int i = 0; i < nEnerBandsUsed; i++){
		for(int j = 0; j < nFramesUsed; j++){
			if(decorMat[j][i] < 0) outputDecorFile << -1 << endl;
			else if(decorMat[j][i] > 0) outputDecorFile << 1 << endl;
			else if(decorMat[j][i] == 0) outputDecorFile << 0 << endl;
			//outputDecorFile << decorMat[j][i] << endl;
		}
	}
	outputDecorFile.close();

	//More Plots

// cout << "ok" << endl;
	/*subplot(3,1,3);
	vector<double> f_axis = {};
	for(unsigned int i = 0; i < (N/2+1)*frameMatrix.size(); i++){
		f_axis.push_back(sampleRate * i / N);		
	}
	plot(f_axis,outV);
	xlim(-0.05,0.05);
	ylim(0,100);
    show();*/
/*
	vector<vector<double>> frameMatrixFft = {};			
	vector<double> fmfFrame = {};
	for(unsigned int i = 0; i < frameMatrix.size(); i++){
		for(int j = 0; j < N/2+1; j++){
			fmfFrame.push_back(outV[i][j]);
		}
		frameMatrixFft.push_back(fmfFrame);
		fmfFrame.clear();
	}
	vector<vector<double>> x,y,z = {};
	
	for(unsigned int i = 0; i < frameMatrix.size(); i+=10){
		vector<double> xrow, yrow, zrow = {};
		for(int j = 0; j < N/2+1; j+=10){
			xrow.push_back(i);
            yrow.push_back(j);
            zrow.push_back(frameMatrixFft[i][j]);
		}
		x.push_back(xrow);
        y.push_back(yrow);
        z.push_back(zrow);
	}*/
	/*for (double i = -5; i <= 5;  i += 0.25) {
        std::vector<double> xrow, yrow, zrow;
        for (double j = -5; j <= 5; j += 0.25) {
            xrow.push_back(i);
            yrow.push_back(j);
            zrow.push_back(::std::sin(::std::hypot(i, j)));
        }
        x.push_back(xrow);
        y.push_back(yrow);
        z.push_back(zrow);
    }*/

	//plot_surface(x,y,z);
	show();


	return 0;
}
