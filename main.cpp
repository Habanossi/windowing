#include "matplotlibcpp.h"
#include <iostream>
#include <stdlib.h>
#include <vector>
#include "string.h"
#include "windowing.cpp"
#include "Wavreader.cpp"
#include <fftw3.h>
#include "resample.h"

using namespace std;

int main(){

	//1.1 Read the audio file and sampling rate
	int sampleRate = readWav();
	int sampleRateTarget = 16000;
	int resampled = 0;
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
	vector<double> dataResampled = {};
	if(sampleRate != sampleRateTarget){
		resample(sampleRateTarget,sampleRate,dataRaw,dataResampled);
		resampled = 1;
	}
	vector<double> dataV = {};
	if(resampled){
		 for(auto i = dataResampled.begin(); i != dataResampled.end(); i++) dataV.push_back(*i);
	}
	else{
		for(auto i = dataRaw.begin(); i != dataRaw.end(); i++) dataV.push_back(*i);
	}


	//1.3 Split the data sequence into windows, implement windowing.cpp
	int frameLength = sampleRate * 0.03;  //25ms in samples, (400?) 			//flattop : 0.03
	int hopSize 	= sampleRate * 0.02; //12.5 ms in samples (50% overlap)	//flattop : 0.02 (2/3 overlap)
	//vector<string> windowTypes = ("rect", "hann", "cosine", "hamming");
	vector<vector<double>> frameMatrix = window("data/raw.dat", frameLength, hopSize, "flattop"); //Windowing

	

	//TODO 1.4 Visualization?
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

	//FFTW
	double matrixSize = frameMatrix.size()*frameLength;					//total size of frameMatrix (n x m)
	double* in = (double*) fftw_malloc(sizeof(double) * matrixSize);	//allocate fft-input array
	fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * matrixSize/2 + 1);	//allocate fft-output array
    fftw_plan p;	//define fft plan
	int N = 1024;
	ofstream offtfile;
	offtfile.open("data/fft.dat"); 	//open textfile for data to matlab for testing
	vector<vector<double>> outV = {};		//vector for fft-output data
	vector<double> outVCol = {};
    p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_MEASURE); //init fft-plan
	for(unsigned int i = 0; i < frameMatrix.size(); i++){		//for each column in frameMatrix, set it as fft-input
		for(int j = 0; j < frameLength; j++){			//
			in[j] = frameMatrix[i][j];					//
		}												//
		fftw_execute(p);								//execute fft for each column
		for(int i = N/2+1; i > 0; i--){												//
			offtfile << 20*log10(abs(out[i][0] + out[i][1])) << endl;								//send output to textfile
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
	int nfft = N/2 + 1;		//513
	int nframes = outV.size(); //216
	int bandsPerEner = floor(nfft/nEnerBands);
	unsigned int count = 0;
	int sum = 0;
	int mean = 0;

	for(int i = nfft-1; i > 0; i--){
		if(count % bandsPerEner == 0){
			for(int j = 0; j < nframes; j++){
				for (int k = 0; k < bandsPerEner; k++){
					if(i-k >= 0) sum += outV[j][i-k]; 
					
				} 
				sum = 10*log10(pow(abs(sum),2));
				mean += sum;										
				enerFrame.push_back(sum);													
				sum = 0;
			}	
			mean /= enerFrame.size();
			for(unsigned int k = 0; k < enerFrame.size(); k++){			
				enerFrame[k] -= mean;
			}	
			enerOut.push_back(enerFrame);
			enerFrame.clear();	
			mean = 0;
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

	//decor_fp run - Extract fingerprint using Mahalanobis decorrelation or linear prediction ( defined by transformation matrix) over context area.

	/*ctx_bands = np.zeros((len(spect_env[0,:]) - x_lim_over - x_lim_under, len(ctx_shape), nener_bands - y_lim_over - y_lim_under))
	
	mult_aux = np.zeros((len(spect_env[0,:]) - x_lim_over - x_lim_under, len(ctx_shape), nener_bands - y_lim_over - y_lim_under))

	decor_mat = np.zeros((len(spect_env[0,:]) - x_lim_over - x_lim_under, nener_bands - y_lim_over - y_lim_under))*/
	ifstream iTMatFile, iCtxOptFile, iCtxShapeFile;
	iTMatFile.open("data/T_mat.txt");
	iCtxOptFile.open("data/ctx_opt.txt");	
	iCtxShapeFile.open("data/ctx_shape.txt");

	vector<double> tMatData, ctxOptData, ctxShapeData = {};
	double a, b, c;
	
	while(iTMatFile >> a){
		tMatData.push_back(a);
	}
	while(iCtxOptFile >> b){
		ctxOptData.push_back(b);
	}
	int d = 0;
	while(iCtxShapeFile >> c){
		ctxShapeData.push_back(c);
		d++;
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
	int nEnerbandsUsed = nEnerBands - yMin - yMax;
	int nFramesUsed = enerOut[0].size() - xMin - xMax;

	vector<vector<vector<double>>> ctxBands(nFramesUsed, vector<vector<double>>(ctxShapeData.size(), vector<double>(nEnerbandsUsed)));
	vector<vector<vector<double>>> multAux (nFramesUsed, vector<vector<double>>(ctxShapeData.size(), vector<double>(nEnerbandsUsed)));
	vector<vector<double>> 		   decorMat(nFramesUsed, vector<double>(nEnerbandsUsed));

	for(int i = 0; i < nEnerBandsUsed; i++){
		for(int j = 0; j < nFramesUsed; j++){
			for(int k = 0; k < ctxShapeData.size(); k++){			
				int ctxIdx = ctxOptData[k*nEnerBandsUsed + j] - 1;
				int idx = j + xMin + ctxShapeData[ctxIdx*2 + 1];
				int idy = i + yMin + ctxShapeData[ctxIdx*2];
				ctxBands[j][k][i] = enerOut[idy][idx]; 
			}
			
		}
	}
	
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
