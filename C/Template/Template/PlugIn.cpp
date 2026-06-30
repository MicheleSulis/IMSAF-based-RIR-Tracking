#include "StdAfx.h"
#include ".\plugin.h"


PlugIn::PlugIn(InterfaceType _CBFunction,void * _PlugRef,HWND ParentDlg): LEEffect(_CBFunction,_PlugRef,ParentDlg)
{
	// Chiamata alla CBFunction
	FrameSize = CBFunction(this,NUTS_GET_FS_SR,0,(LPVOID)AUDIOPROC);
	SampleRate = CBFunction(this,NUTS_GET_FS_SR,1,(LPVOID)AUDIOPROC);	
	
	memset(save_name, 0, MAX_FILE_NAME_LENGTH * sizeof(char));
	delta_h = 0.001;
	mu_h = 0.32;
	P = 2;
	M = 2;
	L = 2;
	D = 2;
	I = 8;
	filter_len = 129;
	K = 128;

	taps = 0;
	h_ana_complex = 0;
	dft_rot_ana = 0;
	dft_rot_syn = 0;
	alpha_rad = 0;
	w_i = 0;
	overlap_x = 0;
	overlap_y = 0;
	H_sub = 0;
	s_state = 0;
	s_head_idx = 0;
	state_syn = 0;
	S_memory = 0;
	tmp_x_full = 0;
	tmp_y_full = 0;
	s_base = 0;
	s_decorr = 0;
	y_sub = 0;
	x_out_frame = 0;
	s_ij = 0;
	u_ij = 0;
	y_pred = 0;
	e_ki = 0;
	S_mat_temp = 0;
	a_i = 0;
	R_mat = 0;
	P_vec = 0;
	overlap_out_x = 0;
	tmp_out_x_full = 0;
	request_rir_save = false;
}

int __stdcall PlugIn::LEPlugin_Process(PinType** Input, PinType** Output, LPVOID ExtraInfo)
{
	int sub_frames = FrameSize / D;

	// TO-DO: implementare i pesi
	// Andrebbero calcolati ad ogni frame, valutare se servono perché le performance sembrano già buone
	if (w_i[0] == 0.0) {
		for (int i = 0; i < I; i++) w_i[i] = 1.0;
	}

	// Massimo 32 speaker e 32 microfoni
	double* in_x[32];
	double* in_y[32];
	double* out_x[32];

	for (int l = 0; l < L; l++) {
		in_x[l] = (double*)Input[l]->DataBuffer;
		out_x[l] = (double*)Output[l]->DataBuffer;
	}
	for (int m = 0; m < M; m++) {
		in_y[m] = (double*)Input[L + m]->DataBuffer;
	}

	// Preparazione buffer overlap per l'ingresso IN_X
	for (int l = 0; l < L; l++) {
		int offset_l = l * (filter_len - 1);
		int tmp_offset = l * (FrameSize + filter_len - 1);

		ippsCopy_64f(&overlap_x[offset_l], &tmp_x_full[tmp_offset], filter_len - 1);
		ippsCopy_64f(in_x[l], &tmp_x_full[tmp_offset + filter_len - 1], FrameSize);
		ippsCopy_64f(&in_x[l][FrameSize - (filter_len - 1)], &overlap_x[offset_l], filter_len - 1);
	}

	// Preparazione buffer overlap per l'ingresso IN_Y
	for (int m = 0; m < M; m++) {
		int offset_m = m * (filter_len - 1);
		int tmp_offset = m * (FrameSize + filter_len - 1);

		ippsCopy_64f(&overlap_y[offset_m], &tmp_y_full[tmp_offset], filter_len - 1);
		ippsCopy_64f(in_y[m], &tmp_y_full[tmp_offset + filter_len - 1], FrameSize);
		ippsCopy_64f(&in_y[m][FrameSize - (filter_len - 1)], &overlap_y[offset_m], filter_len - 1);
	}


	// Analisi di x e decorrelazione
	int temp_global_k = global_k;
	for (int k_local = 0; k_local < sub_frames; k_local++) {
		int n_full = k_local * D;
		int k_mod = temp_global_k % I;

		for (int i = 0; i < I; i++) {
			Ipp64fc rotor_ana = dft_rot_ana[i * I + k_mod];
			for (int l = 0; l < L; l++) {
				Ipp64fc filter_out = { 0.0, 0.0 };
				int tmp_offset = l * (FrameSize + filter_len - 1);

				for (int v = 0; v < filter_len; v++) {
					double sample = tmp_x_full[tmp_offset + n_full + (filter_len - 1) - v];
					Ipp64fc h_val = h_ana_complex[i * filter_len + v];
					filter_out.re += sample * h_val.re;
					filter_out.im += sample * h_val.im;
				}

				Ipp64fc s_base_val;
				ippsMulC_64fc(&filter_out, rotor_ana, &s_base_val, 1);
				s_base[(l * I + i) * sub_frames + k_local] = s_base_val;

				double f_l_val = 2.0 * l;
				double phi = alpha_rad[i] * sin(2.0 * M_PI * f_l_val * (double)temp_global_k * D / SampleRate);
				Ipp64fc phase_mod = { cos(phi), sin(phi) };
				ippsMulC_64fc(&s_base_val, phase_mod, &s_decorr[(l * I + i) * sub_frames + k_local], 1);
			}
		}
		temp_global_k++;
	}

	// Sintesi degli speaker (out_x)
	int temp_global_n = global_n;
	for (int n = 0; n < FrameSize; n++) {
		int n_mod = temp_global_n % I;
		for (int l = 0; l < L; l++) {
			double out_val = 0.0;
			for (int i = 0; i < I; i++) {
				int state_idx = l * I * (filter_len - 1) + i * (filter_len - 1);
				Ipp64fc w_val = { 0.0, 0.0 };

				if (n % D == 0) {
					int k_local = n / D;
					w_val = s_decorr[(l * I + i) * sub_frames + k_local];
				}

				Ipp64fc filter_out = { w_val.re * taps[0], w_val.im * taps[0] };
				for (int v = 1; v < filter_len; v++) {
					Ipp64fc state_val = state_syn[state_idx + (v - 1)];
					filter_out.re += state_val.re * taps[v];
					filter_out.im += state_val.im * taps[v];
				}

				Ipp64fc s_i;
				Ipp64fc rotor_syn = dft_rot_syn[i * I + n_mod];
				ippsMulC_64fc(&filter_out, rotor_syn, &s_i, 1);
				out_val += s_i.re;

				// Shift linea di ritardo
				if (filter_len > 2) ippsMove_64fc(&state_syn[state_idx], &state_syn[state_idx + 1], filter_len - 2);
				if (filter_len > 1) state_syn[state_idx] = w_val;
			}
			out_x[l][n] = out_val * D;
		}
		temp_global_n++;
	}

	// Aggiornamento del filtro adattativo
	// Overlap di out_x
	for (int l = 0; l < L; l++) {
		int offset_l = l * (filter_len - 1);
		int tmp_offset = l * (FrameSize + filter_len - 1);
		ippsCopy_64f(&overlap_out_x[offset_l], &tmp_out_x_full[tmp_offset], filter_len - 1);
		ippsCopy_64f(out_x[l], &tmp_out_x_full[tmp_offset + filter_len - 1], FrameSize);
		ippsCopy_64f(&out_x[l][FrameSize - (filter_len - 1)], &overlap_out_x[offset_l], filter_len - 1);
	}

	for (int k_local = 0; k_local < sub_frames; k_local++) {
		int n_full = k_local * D;
		int k_mod = global_k % I;

		for (int i = 0; i < I; i++) {
			Ipp64fc rotor_ana = dft_rot_ana[i * I + k_mod];

			// Analisi microfono (in_y -> y_sub)
			for (int m = 0; m < M; m++) {
				Ipp64fc filter_out = { 0.0, 0.0 };
				int tmp_offset = m * (FrameSize + filter_len - 1);
				for (int v = 0; v < filter_len; v++) {
					double sample = tmp_y_full[tmp_offset + n_full + (filter_len - 1) - v];
					Ipp64fc h_val = h_ana_complex[i * filter_len + v];
					filter_out.re += sample * h_val.re;
					filter_out.im += sample * h_val.im;
				}
				ippsMulC_64fc(&filter_out, rotor_ana, &y_sub[(m * I + i) * sub_frames + k_local], 1);
			}

			// Analisi speaker reale (out_x -> s_true_val) per usarlo come regressore
			for (int l = 0; l < L; l++) {
				Ipp64fc filter_out = { 0.0, 0.0 };
				int tmp_offset = l * (FrameSize + filter_len - 1);
				for (int v = 0; v < filter_len; v++) {
					double sample = tmp_out_x_full[tmp_offset + n_full + (filter_len - 1) - v];
					Ipp64fc h_val = h_ana_complex[i * filter_len + v];
					filter_out.re += sample * h_val.re;
					filter_out.im += sample * h_val.im;
				}
				Ipp64fc s_true_val;
				ippsMulC_64fc(&filter_out, rotor_ana, &s_true_val, 1);

				// Popola lo stato del regressore con l'analisi di out_x
				int ch_idx = i * L + l;
				s_state[ch_idx * Ki + s_head_idx[ch_idx]] = s_true_val;
				s_head_idx[ch_idx] = (s_head_idx[ch_idx] + 1) % Ki;

				for (int delay = 0; delay < Ki; delay++) {
					int read_idx = (s_head_idx[ch_idx] - 1 - delay + Ki) % Ki;
					s_ij[l * Ki + delay] = s_state[ch_idx * Ki + read_idx];
				}
			}
			double norm_s_val;
			ippsNorm_L2_64fc64f(s_ij, L * Ki, &norm_s_val);
			double norm_s = norm_s_val * norm_s_val;

			// Evita divisioni matematiche per 0 dovute a silenzi in segnali digitali
			if (norm_s < 1e-12) {
				if (P > 0) {
					// Shift della memoria
					ippsMove_64fc(&S_memory[i * P * L * Ki],
						&S_memory[i * P * L * Ki + L * Ki],
						(P - 1) * L * Ki);
					ippsCopy_64fc(s_ij, &S_memory[i * P * L * Ki], L * Ki);
				}
				continue;
			}

			// Decorrelazione
			ippsCopy_64fc(s_ij, u_ij, L * Ki); // Base: u_ij = s_ij

			if (P > 0) {
				double trace_R = 0.0;

				// Costruzione Matrice di Autocorrelazione R_mat e Vettore P_vec
				for (int p1 = 0; p1 < P; p1++) {
					Ipp64fc* S_p1 = &S_memory[i * P * L * Ki + p1 * L * Ki];

					// P_vec[p1] = S_p1^H * s_ij
					Ipp64fc p_val = { 0.0, 0.0 };
					for (int k = 0; k < L * Ki; k++) {
						p_val.re += S_p1[k].re * s_ij[k].re + S_p1[k].im * s_ij[k].im;
						p_val.im += S_p1[k].re * s_ij[k].im - S_p1[k].im * s_ij[k].re;
					}
					P_vec[p1] = p_val;

					for (int p2 = p1; p2 < P; p2++) {
						Ipp64fc* S_p2 = &S_memory[i * P * L * Ki + p2 * L * Ki];

						// R_mat[p1, p2] = S_p1^H * S_p2
						Ipp64fc r_val = { 0.0, 0.0 };
						for (int k = 0; k < L * Ki; k++) {
							r_val.re += S_p1[k].re * S_p2[k].re + S_p1[k].im * S_p2[k].im;
							r_val.im += S_p1[k].re * S_p2[k].im - S_p1[k].im * S_p2[k].re;
						}

						R_mat[p1 * P + p2] = r_val;
						if (p1 != p2) R_mat[p2 * P + p1] = { r_val.re, -r_val.im }; // Hermitiana
					}
					trace_R += R_mat[p1 * P + p1].re;
				}

				// Regolarizzazione
				double reg = 0.1 * trace_R / P + 1e-6;
				for (int p1 = 0; p1 < P; p1++) R_mat[p1 * P + p1].re += reg;

				// Risolutore Lineare (Decomposizione Cholesky In-Place)
				bool chol_ok = true;
				for (int r = 0; r < P; r++) {
					for (int c = 0; c <= r; c++) {
						Ipp64fc sum = R_mat[r * P + c];
						for (int k = 0; k < c; k++) {
							Ipp64fc L_rk = R_mat[r * P + k];
							Ipp64fc L_ck = R_mat[c * P + k];
							sum.re -= (L_rk.re * L_ck.re + L_rk.im * L_ck.im);
							sum.im -= (L_rk.im * L_ck.re - L_rk.re * L_ck.im);
						}
						if (r == c) {
							if (sum.re <= 0) { chol_ok = false; break; }
							R_mat[r * P + c] = { sqrt(sum.re), 0.0 };
						}
						else {
							double inv_Lcc = 1.0 / R_mat[c * P + c].re;
							R_mat[r * P + c] = { sum.re * inv_Lcc, sum.im * inv_Lcc };
						}
					}
					if (!chol_ok) break;
				}

				// Sostituzione Forward & Backward
				if (chol_ok) {
					// Forward: L * y = P_vec (salvato in P_vec)
					for (int r = 0; r < P; r++) {
						Ipp64fc sum = P_vec[r];
						for (int c = 0; c < r; c++) {
							Ipp64fc L_rc = R_mat[r * P + c];
							Ipp64fc y_c = P_vec[c];
							sum.re -= (L_rc.re * y_c.re - L_rc.im * y_c.im);
							sum.im -= (L_rc.re * y_c.im + L_rc.im * y_c.re);
						}
						P_vec[r] = { sum.re / R_mat[r * P + r].re, sum.im / R_mat[r * P + r].re };
					}
					// Backward: L^H * a = y (salvato in a_i)
					for (int r = P - 1; r >= 0; r--) {
						Ipp64fc sum = P_vec[r];
						for (int c = r + 1; c < P; c++) {
							Ipp64fc L_cr = R_mat[c * P + r]; // L[c,r]
							Ipp64fc a_c = a_i[c];
							// sum -= L^H[r,c] * a_c = L[c,r]^* * a_c
							sum.re -= (L_cr.re * a_c.re + L_cr.im * a_c.im); // Coniugato: img invertita
							sum.im -= (L_cr.re * a_c.im - L_cr.im * a_c.re);
						}
						a_i[r] = { sum.re / R_mat[r * P + r].re, sum.im / R_mat[r * P + r].re };
					}

					// Calcolo del segnale Sbiancato: u_ij = s_ij - S_mat * a_i
					for (int p = 0; p < P; p++) {
						Ipp64fc neg_a = { -a_i[p].re, -a_i[p].im };
						Ipp64fc* S_p = &S_memory[i * P * L * Ki + p * L * Ki];
						ippsAddProductC_64fc(S_p, neg_a, u_ij, L * Ki);
					}
				}

				// Aggiornamento memoria dei regressori
				ippsMove_64fc(&S_memory[i * P * L * Ki],
					&S_memory[i * P * L * Ki + L * Ki],
					(P - 1) * L * Ki);
				ippsCopy_64fc(s_ij, &S_memory[i * P * L * Ki], L * Ki);
			}

			// Filtraggio adattativo
			double norm_u_val;
			ippsNorm_L2_64fc64f(u_ij, L * Ki, &norm_u_val);
			double norm_u = norm_u_val * norm_u_val;

			Ipp64fc* H_sub_i = &H_sub[i * M * L * Ki];

			for (int m = 0; m < M; m++) {
				Ipp64fc y_hat = { 0.0, 0.0 };
				Ipp64fc* H_row = &H_sub_i[m * L * Ki];

				// Prodotto scalare senza coniugazione
				ippsDotProd_64fc(H_row, s_ij, L * Ki, &y_hat);

				Ipp64fc y_actual = y_sub[(m * I + i) * sub_frames + k_local];
				Ipp64fc e;
				ippsSub_64fc(&y_hat, &y_actual, &e, 1);

				double step_size = mu_h * w_i[i] / (norm_u + delta_h);
				Ipp64fc update_factor = { e.re * step_size, e.im * step_size };

				// Aggiornamento pesi esplicito
				for (int lk = 0; lk < L * Ki; lk++) {
					// H = H + update_factor * u_ij^*
					H_row[lk].re += update_factor.re * u_ij[lk].re + update_factor.im * u_ij[lk].im;
					H_row[lk].im += update_factor.im * u_ij[lk].re - update_factor.re * u_ij[lk].im;
				}
			}
		} // Fine ciclo Sottobande

		global_k++;
	} // Fine ciclo Decimazione (k_local)

	// Aggiornamento del contatore di sintesi
	global_n = temp_global_n;

	if (request_rir_save) {
		time_t now = time(0);
		struct tm tstruct;
		localtime_s(&tstruct, &now);
		char timestamp[80];
		strftime(timestamp, sizeof(timestamp), "%Y%m%d_%H%M%S", &tstruct);

		sprintf(rir_filename, "%s_EstimatedRIR.dat", timestamp);
		SaveEstimatedRIR(rir_filename);

		request_rir_save = false;
	}

	return COMPLETED;
}

void __stdcall PlugIn::LEPlugin_Init()
{
	Ki = ceil((double)(K + filter_len) / D);

	init_vector(taps, filter_len);
	init_vector(h_ana_complex, I * filter_len);
	init_vector(dft_rot_ana, I * I);
	init_vector(dft_rot_syn, I * I);
	init_vector(alpha_rad, I);
	init_vector(w_i, I);
	global_k = 0;
	global_n = 0;
	init_vector(overlap_x, L * (filter_len - 1));
	init_vector(overlap_y, M * (filter_len - 1));
	init_vector(state_syn, L * I * (filter_len - 1));
	init_vector(H_sub, I * L * M * Ki);
	init_vector(s_state, I * L * Ki);
	if (s_head_idx == 0) {
		s_head_idx = (int*)malloc(I * L * sizeof(int));
		memset(s_head_idx, 0, I * L * sizeof(int));
	}
	init_vector(S_memory, I * P * L * Ki);
	init_vector(tmp_x_full, L*(FrameSize + filter_len - 1));
	init_vector(tmp_y_full, M*(FrameSize + filter_len - 1));
	init_vector(s_base, L * I * FrameSize / D);
	init_vector(s_decorr, L * I * FrameSize / D);
	init_vector(y_sub, M * I * FrameSize / D);
	init_vector(x_out_frame, L * FrameSize);
	init_vector(s_ij, L * Ki);
	init_vector(u_ij, L * Ki);
	init_vector(y_pred, M);
	init_vector(e_ki, M);
	init_vector(S_mat_temp, P * L * Ki);
	init_vector(a_i, P);
	init_vector(R_mat, P * P);
	init_vector(P_vec, P);
	init_vector(overlap_out_x, L * (filter_len - 1));
	init_vector(tmp_out_x_full, L * (FrameSize + filter_len - 1));

	// Caricamento dei tappi del filtro prototipo
	read_dat(save_name, taps, filter_len);

	// Costruzione degli angoli per la modulazione di fase
	for (int i = 0; i < I; i++) {
		double alpha_deg = 0.0;
		if (i <= 3) alpha_deg = 20.0;
		else if (i == 4) alpha_deg = 40.0;
		else if (i == 5) alpha_deg = 70.0;
		else if (i == 6) alpha_deg = 90.0;
		else alpha_deg = 180.0;
		alpha_rad[i] = alpha_deg * (M_PI / 180.0);
	}

	// Calcolo dei filtri complessi modulati
	for (int i = 0; i < I; i++) {
		for (int m = 0; m < filter_len; m++) {
			double angle = (2.0 * M_PI * i * m) / I;

			double real_part = taps[m] * cos(angle);
			double imag_part = taps[m] * sin(angle);

			h_ana_complex[i * filter_len + m] =  Ipp64fc {real_part, imag_part};
		}

		// Calcolo dei coefficienti di fase per il filtro di analisi
		for (int k_mod = 0; k_mod < I; k_mod++) {
			double angle = -(2.0 * M_PI * i * (k_mod * D) / I);
			double real_part = cos(angle);
			double imag_part = sin(angle);
			dft_rot_ana[i * I + k_mod] = Ipp64fc{ real_part, imag_part };
		}

		// Calcolo dei coefficienti di fase per il filtro di sintesi
		for (int n_mod = 0; n_mod < I; n_mod++) {
			double angle = (2.0 * M_PI * i * n_mod) / I;
			double real_part = cos(angle);
			double imag_part = sin(angle);
			dft_rot_syn[i * I + n_mod] = Ipp64fc{ real_part, imag_part };
		}
	}

}

void __stdcall PlugIn::LEPlugin_Delete()
{
	// Salvataggio automatico della RIR stimata prima di eliminare il NUTS
	time_t now = time(0);
	struct tm tstruct;
	localtime_s(&tstruct, &now);
	char timestamp[80];
	strftime(timestamp, sizeof(timestamp), "%Y%m%d_%H%M%S", &tstruct);
	char auto_rir_filename[MAX_FILE_NAME_LENGTH];
	sprintf_s(auto_rir_filename, MAX_FILE_NAME_LENGTH, "%s_Auto_EstimatedRIR.dat", timestamp);
	SaveEstimatedRIR(auto_rir_filename);

	destroy_vector(taps);
	destroy_vector(h_ana_complex);
	destroy_vector(dft_rot_ana);
	destroy_vector(dft_rot_syn);
	destroy_vector(alpha_rad);
	destroy_vector(w_i);
	destroy_vector(overlap_x);
	destroy_vector(overlap_y);
	destroy_vector(state_syn);
	destroy_vector(H_sub);
	if (s_head_idx != 0) {
		free(s_head_idx);
		s_head_idx = 0;
	}
	destroy_vector(S_memory);
	destroy_vector(tmp_x_full);
	destroy_vector(tmp_y_full);
	destroy_vector(s_base);
	destroy_vector(s_decorr);
	destroy_vector(y_sub);
	destroy_vector(x_out_frame);
	destroy_vector(s_ij);
	destroy_vector(u_ij);
	destroy_vector(y_pred);
	destroy_vector(e_ki);
	destroy_vector(S_mat_temp);
	destroy_vector(a_i);
	destroy_vector(P_vec);
	destroy_vector(overlap_out_x);
	destroy_vector(tmp_out_x_full);
}

PlugIn::~PlugIn(void)
{

}


bool __stdcall PlugIn::LEInfoIO(int index,int type, char *StrInfo)
{
	if(type==INPUT && index<L) sprintf(StrInfo,"ToSpeaker[%d]",index);
	if(type==INPUT && ((index>=L) && (index<M+L))) sprintf(StrInfo, "FromMicrophone[%d]", index-L);
	if(type==OUTPUT) sprintf(StrInfo,"Out[%d]",index);
	return true;
}

int __stdcall PlugIn::LESetDefPin(int index,int type, PinType *Info)
{
	if (type==OUTPUT) 
	{
		Info->DataType=PLAYBUFFER;
		Info->Exclusive=true; 
		Info->DataLen=FrameSize;
		Info->MaxDataLen=FrameSize;
		return OUTPUT; 
	}

	if (type==INPUT) 
	{
		Info->DataType=PLAYBUFFER;
		Info->Exclusive=true; 
		Info->DataLen=FrameSize;
		Info->MaxDataLen=FrameSize;
		return INPUT; 
	}

	return -1;
}

int  __stdcall PlugIn::LEConnectionRequest(int IOType,int Index,PinType *NewType)
{
	return 0;
}

LPVOID __stdcall PlugIn::LEOnNUTechMessage(int MessageType,int MessageID,WPARAM wParam,LPARAM lParam)
{
	return 0;
}


void __stdcall PlugIn::LESetName(char *Name)
{
	strncpy(Name, NUTS_NAME, MAXNAME);
}

void __stdcall PlugIn::LESetParameter(int Index,void *Data,LPVOID bBroadCastInfo)
{
	if (Index == IDX_FILTER_LEN) {
		int new_filter_len = *(int*)Data;
		if (new_filter_len != filter_len) {
			filter_len = new_filter_len;
			CBFunction(this, NUTS_UPDATERTWATCH, IDX_FILTER_LEN, 0);
		}
		return;
	}

	if (Index == IDX_SAVE_NAME) {
		char* new_save_name = (char*)Data;
		if (strcmp(new_save_name, save_name) != 0) {
			strncpy(save_name, new_save_name, MAX_FILE_NAME_LENGTH * sizeof(char));
			CBFunction(this, NUTS_UPDATERTWATCH, IDX_SAVE_NAME, 0);
		}
		return;
	}

	if (Index == IDX_L) {
		int new_L = *(int*)Data;
		if ((new_L != L) && (new_L <= 32)) {
			L = new_L;
			LESetNumInput(L + M);
			CBFunction(this, NUTS_SETNUMINPUT, L + M, 0);
			LESetNumOutput(L);
			CBFunction(this, NUTS_SETNUMOUTPUT, L, 0);
			CBFunction(this, NUTS_UPDATERTWATCH, IDX_L, 0);
		}
		return;
	}

	if (Index == IDX_M) {
		int new_M = *(int*)Data;
		if ((new_M != M) && (new_M <= 32)) {
			M = new_M;
			LESetNumInput(L + M);
			CBFunction(this, NUTS_SETNUMINPUT, L + M, 0);
			CBFunction(this, NUTS_UPDATERTWATCH, IDX_M, 0);
		}
		return;
	}

	if (Index == IDX_P) {
		int new_P = *(int*)Data;
		if (new_P != P) {
			P = new_P;
			CBFunction(this, NUTS_UPDATERTWATCH, IDX_P, 0);
		}
		return;
	}

	if (Index == IDX_DELTA_H) {
		double new_delta_h = *(double*)Data;
		if (new_delta_h != delta_h) {
			delta_h = new_delta_h;
			CBFunction(this, NUTS_UPDATERTWATCH, IDX_DELTA_H, 0);
		}
		return;
	}

	if (Index == IDX_MU_H) {
		double new_mu_h = *(double*)Data;
		if (new_mu_h != mu_h) {
			mu_h = new_mu_h;
			CBFunction(this, NUTS_UPDATERTWATCH, IDX_MU_H, 0);
		}
		return;
	}

	if (Index == IDX_SUBBAND_NUM) {
		int new_I = *(int*)Data;
		if (new_I != I) {
			I = new_I;
			CBFunction(this, NUTS_UPDATERTWATCH, IDX_SUBBAND_NUM, 0);
		}
		return;
	}

	if (Index == IDX_DEC_FACTOR) {
		int new_D = *(int*)Data;
		if (new_D != D) {
			D = new_D;
			CBFunction(this, NUTS_UPDATERTWATCH, IDX_DEC_FACTOR, 0);
		}
		return;
	}

	if (Index == IDX_RIR_SAMPLES) {
		int new_samples = *(int*)Data;
		if (new_samples != K) {
			K = new_samples;
			CBFunction(this, NUTS_UPDATERTWATCH, IDX_RIR_SAMPLES, 0);
		}
		return;
	}

	if (Index == IDX_SAVE_RIR_TRIGGER) {
		int trigger = *(int*)Data;
		if (trigger == 1) {
			request_rir_save = true;
			int reset_val = 0;
			CBFunction(this, NUTS_UPDATERTWATCH, IDX_SAVE_RIR_TRIGGER, &reset_val);
		}
		return;
	}

}

int  __stdcall PlugIn::LEGetParameter(int Index,void *Data)
{
	if (Index == IDX_FILTER_LEN) {
		*(int*)Data = filter_len;
		return 0;
	}

	if (Index == IDX_SAVE_NAME) {
		strncpy((char*)Data, save_name, MAX_FILE_NAME_LENGTH * sizeof(char));
		return 0;
	}

	if (Index == IDX_L) {
		*(int*)Data = L;
		return 0;
	}

	if (Index == IDX_M) {
		*(int*)Data = M;
		return 0;
	}

	if (Index == IDX_SUBBAND_NUM) {
		*(int*)Data = I;
		return 0;
	}

	if (Index == IDX_DEC_FACTOR) {
		*(int*)Data = D;
		return 0;
	}

	if (Index == IDX_DELTA_H) {
		*(double*)Data = delta_h;
		return 0;
	}

	if (Index == IDX_P) {
		*(int*)Data = P;
		return 0;
	}

	if (Index == IDX_MU_H) {
		*(double*)Data = mu_h;
		return 0;
	}

	if (Index == IDX_RIR_SAMPLES) {
		*(int*)Data = K;
		return 0;
	}


	return 0;
}

void __stdcall PlugIn::LESaveSetUp()
{
	
}

void __stdcall PlugIn::LELoadSetUp()
{

}

void __stdcall PlugIn::LERTWatchInit()
{
	// Watch per la lunghezza del filtro
	WatchType FilWatch;
	memset(&FilWatch, 0, sizeof(WatchType));
	FilWatch.EnableWrite = true;
	FilWatch.LenByte = sizeof(int);
	FilWatch.TypeVar = WTC_INT;
	FilWatch.IDVar = IDX_FILTER_LEN;
	sprintf(FilWatch.VarName, "Filter Length");
	CBFunction(this, NUTS_ADDRTWATCH, 0, (LPVOID)&FilWatch);

	// Watch per il nome del file dei tappi
	WatchType NameWatch;
	memset(&NameWatch, 0, sizeof(WatchType));
	NameWatch.EnableWrite = true;
	NameWatch.LenByte = MAX_FILE_NAME_LENGTH * sizeof(char);
	NameWatch.TypeVar = WTC_LPCHAR;
	NameWatch.IDVar = IDX_SAVE_NAME;
	sprintf_s(NameWatch.VarName, MAXCARDEBUGPLUGIN, "Filter directory");
	CBFunction(this, NUTS_ADDRTWATCH, 0, (LPVOID)&NameWatch);

	// Watch per il numero di altoparlanti
	WatchType LWatch;
	memset(&LWatch, 0, sizeof(WatchType));
	LWatch.EnableWrite = true;
	LWatch.LenByte = sizeof(int);
	LWatch.TypeVar = WTC_INT;
	LWatch.IDVar = IDX_L;
	sprintf(LWatch.VarName, "Number of speakers L");
	CBFunction(this, NUTS_ADDRTWATCH, 0, (LPVOID)&LWatch);

	// Watch per il numero di microfoni
	WatchType MWatch;
	memset(&MWatch, 0, sizeof(WatchType));
	MWatch.EnableWrite = true;
	MWatch.LenByte = sizeof(int);
	MWatch.TypeVar = WTC_INT;
	MWatch.IDVar = IDX_M;
	sprintf(MWatch.VarName, "Number of microphones M");
	CBFunction(this, NUTS_ADDRTWATCH, 0, (LPVOID)&MWatch);

	// Watch per il numero di campioni della rir stimata
	WatchType RIRWatch;
	memset(&RIRWatch, 0, sizeof(WatchType));
	RIRWatch.EnableWrite = true;
	RIRWatch.LenByte = sizeof(int);
	RIRWatch.TypeVar = WTC_INT;
	RIRWatch.IDVar = IDX_RIR_SAMPLES;
	sprintf(RIRWatch.VarName, "Number of samples of the estimated RIR");
	CBFunction(this, NUTS_ADDRTWATCH, 0, (LPVOID)&RIRWatch);

	// Watch per l'ordine di predizione P
	WatchType PWatch;
	memset(&PWatch, 0, sizeof(WatchType));
	PWatch.EnableWrite = true;
	PWatch.LenByte = sizeof(int);
	PWatch.TypeVar = WTC_INT;
	PWatch.IDVar = IDX_P;
	sprintf(PWatch.VarName, "Decorrelation order P");
	CBFunction(this, NUTS_ADDRTWATCH, 0, (LPVOID)&PWatch);

	// Watch per I
	WatchType IWatch;
	memset(&IWatch, 0, sizeof(WatchType));
	IWatch.EnableWrite = true;
	IWatch.LenByte = sizeof(int);
	IWatch.TypeVar = WTC_INT;
	IWatch.IDVar = IDX_SUBBAND_NUM;
	sprintf(IWatch.VarName, "Number of Subbands I");
	CBFunction(this, NUTS_ADDRTWATCH, 0, (LPVOID)&IWatch);

	// Watch per D
	WatchType DWatch;
	memset(&DWatch, 0, sizeof(WatchType));
	DWatch.EnableWrite = true;
	DWatch.LenByte = sizeof(int);
	DWatch.TypeVar = WTC_INT;
	DWatch.IDVar = IDX_DEC_FACTOR;
	sprintf(DWatch.VarName, "Decimation Factor D");
	CBFunction(this, NUTS_ADDRTWATCH, 0, (LPVOID)&DWatch);

	// Watch per lo step size mu_h
	WatchType muWatch;
	memset(&muWatch, 0, sizeof(WatchType));
	muWatch.EnableWrite = true;
	muWatch.LenByte = sizeof(double);
	muWatch.TypeVar = WTC_DOUBLE;
	muWatch.IDVar = IDX_MU_H;
	sprintf(muWatch.VarName, "Step size mu_h");
	CBFunction(this, NUTS_ADDRTWATCH, 0, (LPVOID)&muWatch);

	// Watch per il parametro di regolarizzazione delta_h
	WatchType delta_hWatch;
	memset(&delta_hWatch, 0, sizeof(WatchType));
	delta_hWatch.EnableWrite = true;
	delta_hWatch.LenByte = sizeof(double);
	delta_hWatch.TypeVar = WTC_DOUBLE;
	delta_hWatch.IDVar = IDX_MU_H;
	sprintf(delta_hWatch.VarName, "Step size mu_h");
	CBFunction(this, NUTS_ADDRTWATCH, 0, (LPVOID)&delta_hWatch);

	// Watch per il salvataggio manuale della rir
	WatchType SaveWatch;
	memset(&SaveWatch, 0, sizeof(WatchType));
	SaveWatch.EnableWrite = true;
	SaveWatch.LenByte = sizeof(int);
	SaveWatch.TypeVar = WTC_INT; // (0 = no, 1 = salva)
	SaveWatch.IDVar = IDX_SAVE_RIR_TRIGGER;
	sprintf(SaveWatch.VarName, "TRIGGER: Save RIR to File");
	CBFunction(this, NUTS_ADDRTWATCH, 0, (LPVOID)&SaveWatch);
}

void __stdcall PlugIn::LESampleRateChange(int NewVal,int TrigType)
{
	if(TrigType==AUDIOPROC)
	{
		if(NewVal!=SampleRate)
		{
			SampleRate = NewVal;
		}
	}

} 

void __stdcall PlugIn::LEFrameSizeChange (int NewVal,int TrigType)
{
	if(TrigType==AUDIOPROC)
	{
		if(NewVal!=FrameSize)
		{
			FrameSize = NewVal;
		}
	}
} 

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////     LOADER      ///////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

LEEffect * __stdcall LoadEffect(InterfaceType _CBFunction,void *PlugRef,HWND ParentDlg)
{
	return new PlugIn(_CBFunction,PlugRef,ParentDlg);
}

int __stdcall UnLoadEffect(PlugIn *effect)
{
	delete effect;
	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////     GetStartUpInfo      //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

void __stdcall LENUTSDefProps(char *NameEffect,int *Width, void *data)
{
	
	strncpy(NameEffect,NUTS_NAME,MAXNAME);
	*Width=WIDTHDEF;

	if (data!=0)
	{
		StartUpNUTSProps *Info=(StartUpNUTSProps *)data;
		Info->NumInStartUp=4;
		Info->NumOutStartUp=2;
		Info->BitMaskProc = AUDIOPROC;
		Info->BitMaskDriver = OFFLINEDRIVER | DIRECTDRIVER | ASIODRIVER;
	}
}

void PlugIn::SaveEstimatedRIR(const char* filename)
{
	int impulse_len = K + 2 * filter_len;
	int sub_frames_test = (impulse_len - 1) / D + 1;
	int recon_len = (sub_frames_test - 1) * D + 1;

	Ipp64f* x_test = 0;
	init_vector(x_test, impulse_len);
	x_test[0] = 1.0;

	Ipp64fc* s_sub = 0;
	init_vector(s_sub, I * sub_frames_test);

	Ipp64f* sim_ana_state = 0;
	init_vector(sim_ana_state, filter_len - 1);

	Ipp64fc* y_sub_sim = 0;
	init_vector(y_sub_sim, I * sub_frames_test);

	Ipp64fc* sim_syn_state = 0;
	init_vector(sim_syn_state, I * (filter_len - 1));

	Ipp64f* H_full_path = 0;
	init_vector(H_full_path, recon_len);

	for (int k = 0; k < sub_frames_test; k++) {
		int n_full = k * D;
		for (int i = 0; i < I; i++) {
			Ipp64fc filter_out = { 0.0, 0.0 };
			for (int v = 0; v < filter_len; v++) {
				int read_idx = n_full + (filter_len - 1) - v;
				double sample = 0.0;
				if (read_idx < filter_len - 1) sample = sim_ana_state[read_idx];
				else if (read_idx - (filter_len - 1) < impulse_len) sample = x_test[read_idx - (filter_len - 1)];

				Ipp64fc h_val = h_ana_complex[i * filter_len + v];
				filter_out.re += sample * h_val.re;
				filter_out.im += sample * h_val.im;
			}
			double angle = -(2.0 * M_PI * i * n_full) / I;
			Ipp64fc rotor = { cos(angle), sin(angle) };
			ippsMulC_64fc(&filter_out, rotor, &s_sub[i * sub_frames_test + k], 1);
		}
	}

	std::ofstream outfile(filename, std::ios::binary);
	if (!outfile.is_open()) return;

	int delay_calib = filter_len - 1;

	for (int m = 0; m < M; m++) {
		for (int l = 0; l < L; l++) {

			ippsZero_64fc(y_sub_sim, I * sub_frames_test);
			ippsZero_64fc(sim_syn_state, I * (filter_len - 1));
			ippsZero_64f(H_full_path, recon_len);

			for (int i = 0; i < I; i++) {
				Ipp64fc* H_path = &H_sub[i * M * L * Ki + m * L * Ki + l * Ki];

				for (int k = 0; k < sub_frames_test; k++) {
					Ipp64fc y_val = { 0.0, 0.0 };
					for (int delay = 0; delay < Ki; delay++) {
						if (k - delay >= 0) {
							Ipp64fc s_val = s_sub[i * sub_frames_test + (k - delay)];
							Ipp64fc h_w = H_path[delay];
							y_val.re += h_w.re * s_val.re - h_w.im * s_val.im;
							y_val.im += h_w.re * s_val.im + h_w.im * s_val.re;
						}
					}
					y_sub_sim[i * sub_frames_test + k] = y_val;
				}
			}

			for (int n = 0; n < recon_len; n++) {
				double out_val = 0.0;
				for (int i = 0; i < I; i++) {
					int state_idx = i * (filter_len - 1);

					Ipp64fc w_val = { 0.0, 0.0 };
					if (n % D == 0) w_val = y_sub_sim[i * sub_frames_test + (n / D)];

					Ipp64fc filter_out = { w_val.re * taps[0], w_val.im * taps[0] };
					for (int v = 1; v < filter_len; v++) {
						Ipp64fc state_val = sim_syn_state[state_idx + (v - 1)];
						filter_out.re += state_val.re * taps[v];
						filter_out.im += state_val.im * taps[v];
					}

					double angle = (2.0 * M_PI * i * n) / I;
					Ipp64fc rotor = { cos(angle), sin(angle) };
					out_val += (filter_out.re * rotor.re - filter_out.im * rotor.im);

					if (filter_len > 2) {
						ippsMove_64fc(&sim_syn_state[state_idx], &sim_syn_state[state_idx + 1], filter_len - 2);
					}
					if (filter_len > 1) sim_syn_state[state_idx] = w_val;
				}
				H_full_path[n] = out_val * D;
			}

			for (int k_idx = 0; k_idx < K; k_idx++) {
				double val = 0.0;
				if (delay_calib + k_idx < recon_len) {
					val = H_full_path[delay_calib + k_idx];
				}
				outfile.write(reinterpret_cast<const char*>(&val), sizeof(double));
			}
		}
	}

	outfile.close();
	destroy_vector(x_test);
	destroy_vector(s_sub);
	destroy_vector(sim_ana_state);
	destroy_vector(y_sub_sim);
	destroy_vector(sim_syn_state);
	destroy_vector(H_full_path);
}