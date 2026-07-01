#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include "LEEffect.h"
#include "ipp.h"
#include <numbers>
#include "myLib.h"
#include<fstream>
#include<ctime>

#define PLAYBUFFER_TXT "PlayBuffer"
#define	SCALAR_TXT "Scalar"
#define VECTOR_TXT "Vector"
#define VARLEN_TXT "VarLen"

#define IDX_FILTER_LEN 0
#define IDX_SAVE_NAME 1
#define IDX_NUM_INPUT 2
#define IDX_NUM_OUTPUT 3
#define IDX_SUBBAND_NUM 4
#define IDX_DEC_FACTOR 5
#define IDX_L 6
#define IDX_M 7
#define IDX_DELTA_H 8
#define IDX_P 9
#define IDX_MU_H 10
#define IDX_RIR_SAMPLES 11
#define IDX_SAVE_RIR_TRIGGER 12
			
#define INT16_TXT "Integer 16bit"
#define INT32_TXT "Integer 32bit"
#define FLT32_TXT "Float 32bit"
#define DBL64_TXT "Double 64bit"

#define DATATYPEBITMASK VECTOR|SCALAR|BUFFSTREAMING|VARLEN
#define DATALENBITMASK DATAINT32|DATAINT16|DATAINT8|DATAFLOAT32|DATADOUBLE64|CUSTOM

#define WIDTHDEF 75

#define NUTS_NAME	"IMSAF"
#define MAX_FILE_NAME_LENGTH 255

class PlugIn :	public LEEffect
{
public:
	PlugIn(InterfaceType _CBFunction,void * _PlugRef,HWND ParentDlg);

	~PlugIn(void);
		
	void __stdcall LESetName(char *Name);

	void __stdcall LEPlugin_Init();
	int  __stdcall LEPlugin_Process(PinType **Input,PinType **Output,LPVOID ExtraInfo);
	void __stdcall LEPlugin_Delete(void);

	bool __stdcall LEInfoIO(int index,int type, char *StrInfo);

	int __stdcall LEGetNumInput(){return LEEffect::LEGetNumInput();};  
	int __stdcall LEGetNumOutput(){return LEEffect::LEGetNumOutput();};
	int __stdcall LESetNumInput(int Val,PinType *TypeNewIn=0){return LEGetNumInput();}; 
	int __stdcall LESetNumOutput(int Val,PinType *TypeNewOut=0){return LEGetNumOutput();};

	void __stdcall LESetParameter(int Index,void *Data,LPVOID bBroadCastInfo);
	int  __stdcall LEGetParameter(int Index,void *Data);

	int __stdcall LESetDefPin(int index,int type, PinType *Info);

	HWND __stdcall LEGetWndSet(){return 0;};
	int __stdcall  LEWinSetStatusChange(int NewStatus){return 0;};		

	void __stdcall LESaveSetUp();
	void __stdcall LELoadSetUp();

	void __stdcall LESampleRateChange(int NewVal,int TrigType); 
	void __stdcall LEFrameSizeChange (int NewVal,int TrigType); 
	void __stdcall LEConnectionChange(int IOType,int Index,bool Connected){};
	int  __stdcall LEConnectionRequest(int IOType,int Index,PinType *NewType);
	int  __stdcall LEExtraInfoPinChange(int IOType,int Index,PinExtraInfoType ExInfo){return 0;};

	void __stdcall LERTWatchInit();
	LPVOID __stdcall LEOnNUTechMessage(int MessageType,int MessageID,WPARAM wParam,LPARAM lParam);

private:
	// Nel codice matlab filter_len = V
	int FrameSize,SampleRate, filter_len, L, M, I, D, P, K, Ki;
	// int num_input, num_output;
	char save_name[MAX_FILE_NAME_LENGTH];
	double delta_h, mu_h, fs;
	// Tappi del filtro prototipo
	Ipp64f* taps;
	// Filtri complessi modulati per l'analisi, di lunghezza I*V
	Ipp64fc* h_ana_complex;
	// Parte reale del filtro di analisi invertita
	Ipp64f* h_ana_re_rev;
	// Parte immaginaria del filtro di analisi invertita
	Ipp64f* h_ana_im_rev;
	// Buffer per il coniugato di u_ij
	Ipp64fc* u_ij_conj;
	// esponenziali DFT per l'analisi, di dimensione I*I
	Ipp64fc* dft_rot_ana;
	// esponenziali iDFT per la sintesi, di dimensione I*I
	Ipp64fc* dft_rot_syn;
	// Valori per la decorrelaizone di fase espressi in radianti, di lunghezza I 
	Ipp64f* alpha_rad;
	// Pesi di normalizzazione per le sottobande, di lunghezza I
	Ipp64f* w_i;
	// Indice dei campioni in sottobanda (deve sopravvivere da un frame all'altro)
	unsigned long int global_k;
	// Indice globale dei campioni fullband
	unsigned long int global_n;
	// Memoria per i filtri di analisi che contiene gli ultimi V-1 campioni di input del frame precedente, di dimensione L*(V-1)
	Ipp64f* overlap_x;
	// Memoria per i filtri di analisi dei microfoni, di dimensione M*(V-1)
	Ipp64f* overlap_y;
	// Coda del segnale sovracampionato in sottobanda prima di applicare il filtro prototipo per la sintesi, di lunghezza L*I*(V-1)
	Ipp64fc* state_syn;
	// RIR stimata in sottobanda, di dimensione I*M*L*K_i
	Ipp64fc* H_sub;
	// Linea di ritardo per la regressione, di dimensione I*L*Ki
	Ipp64fc* s_state;
	// Indici per usare s_state come memoria FIFO circolare, di dimensione I*L
	int* s_head_idx;
	// Buffer circolare per la sintesi degli speaker
	int* syn_head;
	// Memoria dei predittori per decorrelazione, di dimensione I*P*L*Ki
	Ipp64fc* S_memory; 

	// Buffer per unire overlap_x e il nuovo frame in ingresso, di dimensione FrameSize+V-1
	Ipp64f* tmp_x_full;
	// Buffer per unire overlap_y con il nuovo ingresso dai microfoni, di dimensione FrameSize+V-1
	Ipp64f* tmp_y_full;
	// Output del filtro di analisi, di dimensione L*I*FrameSize/D
	Ipp64fc* s_base;
	// Segnale in sottobanda dopo la modulazione di fase, di dimensione L*I*FrameSize/D
	Ipp64fc* s_decorr; 
	// Output di analisi dei microfoni, di dimensione M*I*FrameSize/D
	Ipp64fc* y_sub;
	// Sintesi di s_decorr, da scrivere sui pin di uscita per gli altoparlanti, di dimensione L*FrameSize
	Ipp64f* x_out_frame;
	// Regressore combinato all'indice k corrente, di dimensione L*Ki
	Ipp64fc* s_ij;
	// Regressore sbiancato, di dimensione L*Ki
	Ipp64fc* u_ij;
	// Predizione attuale del microfono, di dimensione M
	Ipp64fc* y_pred;
	// Errore per la sottobanda I-esima, di dimensione M
	Ipp64fc* e_ki;

	// Copia di estrazione da S_memory per calcoli, di dimensione P*L*K_i
	Ipp64fc* S_mat_temp;
	// Vettore dei coefficienti del predittore lineare, di dimensione P
	Ipp64fc* a_i;
	// matrice di autocorrelazione S^H * S, di dimensione PxP
	Ipp64fc* R_mat;
	// Vettore di correlazione incrociata S^H * s, di dimensione P
	Ipp64fc* P_vec;
	// Flag per salvare la RIR stimata quando richiesto

	// Buffer per i vettori coniugati in modo da poter usare le ipp per i calcoli
	Ipp64fc* s_ij_conj;
	Ipp64fc* S_memory_conj;

	Ipp64f* overlap_out_x;
	Ipp64f* tmp_out_x_full;
	bool request_rir_save;
	char rir_filename[MAX_PATH];

	void SaveEstimatedRIR(const char* filename);
};