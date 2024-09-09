#include "TaskData.h"

namespace TDVesselWrite {
	using namespace TaskData;
	using namespace Globals;

	void TDWrite(Zadacha& Z);

	inline int WriteCalcData ( Zadacha& Z) {

	    int IsDataWrite, n_wr;
	    double perc;

	    vector<double> out(3);
	    double P;

	    //cout << "T = " << T << " , Time_calc = " << time_calc << " , dt = " << dt << endl;
	    //if ( (T >= write_time)&&((last_write < 0.0000000001 )||(T - last_write > write_delta)))
        if (T - last_write > write_delta)
            IsDataWrite = 1;
        else
            IsDataWrite = 0;

        if (IsDataWrite == 1){

            n_wr = int((T - write_time)/write_delta);
            last_write = n_wr*write_delta + write_time;
            perc = T/time_calc*100;
            cout << "T = " << T << " , Time_calc = " << time_calc << " , dt = " << dt << endl;
            cout << "Progress " << perc << " % " << endl;

            //steps_skipped = 0;
            Globals::Nz++;

            TDWrite(Z); // writing time is inside
            //

            //ofstream fou;
            //fou.open(Z.Pfile, std::ios_base::app);
            //out = TreeLst[0].B[0].URSOB(TreeLst[0].B[0].VB[0][5],TreeLst[0].B[0].VB[1][5]);
            //fou << out[0]/1333.22 << endl;
            //fou.close();

            //fou.open(Z.Qfile, std::ios_base::app);
            //fou << TreeLst[0].B[0].VB[0][5]*TreeLst[0].B[0].VB[1][5] << endl;
            //fou.close();
        }


		return IsDataWrite;
	}
	void TDWriteCurrentTime(string filename) {
		ofstream fout;
		fout.open (filename, ofstream::app);
		fout << Globals::T << endl;
		fout.close();
	}
	void TDWriteCurrentTimeSteps(string filename) {
		ofstream fout;
		fout.open (filename, ofstream::out);
		fout << Globals::Nz << endl;
		fout.close();
	}
//------- ������� ������ � ���� �������-��������� ���������� -------------------------
	void WriteResultLymph ( Derevo& Tr , Zadacha& Z ) {
		double q;
		ofstream fout;
		for ( int i = 0; i < Tr.Nbr; ++i ) {
			fout.open (Tr.B[i].TD.fnameTD[0], ofstream::app);
			for (int k = 0; k < Tr.B[i].pts; ++k) {
				fout << Tr.B[i].VB[0][k] * Tr.B[i].VB[1][k] << " ";
			}
			fout << endl;
			fout.close();
		}
	}
	void WriteResultBronch ( Derevo& Tr , Zadacha& Z ) { }
	void WriteResultArtery ( Derevo& Tr , Zadacha& Z ) { }
	void WriteResultVein ( Derevo& Tr , Zadacha& Z ) { }

	void WriteResultTD ( Derevo& Tr , Zadacha& Z ) {
		if ( Tr.ID == LYMPHATIC ) WriteResultLymph(Tr, Z);
		if ( ( Tr.ID == PULMART ) || ( Tr.ID == SYSART ) ) WriteResultArtery(Tr, Z);
		if ( ( Tr.ID == PULMVEN ) || ( Tr.ID == SYSVEN ) ) WriteResultVein(Tr, Z);
		if ( Tr.ID == BRONCHIAL ) WriteResultBronch(Tr, Z);
	}

//-------������������ ������ � ���� ���������� s,u,p � ������------------------------------
// ����� �������� �����, ������ ��� ��������� ��������
	void WriteResultV(Derevo& Tr , Zadacha& Z) {
		double p;
		ofstream fout;
		vector<double> out(3);

		//int cur_track = 0; // ���� ������� �� ��� �����
		for ( int i = 0; i < Tr.Nbr; i++ ) {
			/*if (Z.N_towrite) {
				if ( i == Z.Num_towrite[cur_track] ) {
					cur_track++;
				} else {
					continue;
				}
			}*/
			for ( int j = 0; j < Z.Cor; ++j ) { 				// s, u
				fout.open (Tr.B[i].TD.fnameVar[j], ofstream::app);
				for (int k = 0; k < Tr.B[i].pts; ++k) {
					fout << Tr.B[i].VB[j][k] << endl;
				}
				fout.close();
			}


			fout.open (Tr.B[i].TD.fnameVar[Z.Cor], ofstream::app); // Z.Cor = 2
			for (int k = 0; k < Tr.B[i].pts; k++) { 			// p
				out =  Tr.B[i].URSOB(Tr.B[i].VB[0][k],Tr.B[i].VB[1][k]);
				p = out[0]/1333.2;
				fout << p << endl;
			}
			fout.close();

			fout.open (Tr.B[i].TD.fnameVar[Z.Cor + 1], ofstream::app);
			for (int k = 0; k < Tr.B[i].pts; k++) { 			// q
				p = Tr.B[i].VB[0][k]*Tr.B[i].VB[1][k];
				fout << p << endl;
			}
			fout.close();

			fout.open (Tr.B[i].TD.fnameVar[Z.Cor + 2], ofstream::app);
			for (int k = 0; k < Tr.B[i].pts; k++) { 			// pave
				p = Tr.B[i].TD.Pave_next/1333.2;
				fout << p << endl;
			}
			fout.close();

        }
		//WriteResultTD(Tr , Z);
		TDGlobals::isFirstTime = 0;
	}

	void WriteBranchResultForGNUPLOT (Derevo& Tr , Zadacha& Z) {}

	void TDWrite(Zadacha& Z) {
	//-------------------------- MATLAB output ---------------------------------
		if ( Globals::use_MATLAB_out ) {
			for ( int i = 0; i < Z.Ntr; ++i ) {
				WriteResultV(TreeLst[i], Z);
			}
			//if ( Z.useLungs ) {
			//	WriteResultV(BronchTree, Z);
			//}
			//cout << "T = " << T << " , Time_calc = " << time_calc << " , dt = " << dt << endl;
			TDWriteCurrentTimeSteps(trim(TDGlobals::root_out) + "tsteps.tres");
			TDWriteCurrentTime(trim(TDGlobals::root_out) + "time.tres");
		}
	//-------------------------- GNUPLOT output ---------------------------------
		/*if ( Globals::use_GNUPLOT_out ) {
			for ( int i = 0; i < Z.Ntr; ++i ) {
				WriteBranchResultForGNUPLOT (TreeLst[i] , Z );
			}
			if ( Z.useLungs ) {
				WriteBranchResultForGNUPLOT ( BronchTree, Z );
			}
		}*/
	}
//===========================================================================================
//  ���������� ������� �������� ��� ������������� � ���������� �����������
//===========================================================================================
	void SaveDump(Derevo& Tr , Zadacha& Z) {
		string fnDump;
		int LResult;

		LResult = MAKEDIRQQ( trim( SharedDirectory ) + trim(Tr.dirname) );
		LResult = MAKEDIRQQ( trim( SharedDirectory ) + trim(Tr.dirname) + slash + "result" );
		fnDump = trim(SharedDirectory) + trim(Tr.dirname) + slash + "result" + slash + "tree" + trim( Achar( Ichar( "0" ) + Tr.ID ) )  +  ".dmp";
		cout << "Dump file : " << fnDump << endl;
		ofstream fout( fnDump , ofstream::out );

		fout <<  Tr.ID << endl;
		for ( int i = 0; i < Tr.Nbr; ++i ) {
			fout <<  Tr.B[ i ].ID << " " <<  Tr.B[ i ].pts << endl;;		// ����� �����, ���������� ����� ���������
			for ( int j = 0; j < Z.Cor; ++j ) {
				for (int k = 0; k < Tr.B[i].pts; ++k) {
					fout << Tr.B[i].VB[j][k] << " ";
				}
				fout << endl;
			}
			for ( int j = 0; j < Z.Cor; ++j ) {
				for (int k = 0; k < Tr.B[i].pts; ++k) {
					fout << Tr.B[i].VBO[j][k] << " ";
				}
				fout << endl;
			}
		}
		fout.close();
	}

	void ConvertResults(Derevo& Tr , Zadacha& Z) {}

	void WriteFFR(Derevo& Tr , Zadacha& Z) {

	    double Pd,Pa, FFR, placeholder;
	    long id_next, addBr, nd;

	    ofstream fou;
        fou.open(Tr.FFRfile);

        placeholder = - 1.5;
        Pa = Tr.B[0].Pave[5];

        for (long i = 0; i < Z.Nsten; i++){

            if (Tr.B[Z.IDsten[i] - 1].stenType == 1){

                id_next = (*(*Tr.B[Z.IDsten[i] - 1].Kn2).Bou[0]).ID;

                Pd = Tr.B[id_next - 1].Pave[2]; // or Pave[0]
                 // near bifurcation
                FFR = Pd/Pa;

                if (id_next <= (Tr.NbrL + 2))
                    addBr = 2;
                else
                    addBr = Tr.NbrL + 2;

                fou <<  Tr.B[Z.IDsten[i] - 1].ID << endl;
                fou <<  FFR << endl;
                fou << Tr.B[Z.IDsten[i] - 1].Qave << endl;
                fou << placeholder << endl;
                //fou << endl;

            }

        }
        cout << "FFR calculated" << endl;

        fou.close();


        fou.open(Tr.FFRfullfile);

        fou <<  "Aorta, mmHg" << endl;
        fou <<  Pa/1333.2 << endl;
        fou <<  "Aorta, Qave" << endl;
        fou <<  Tr.B[0].Qave  << endl;


        fou <<  "LCA" << endl;

        for (long i = 0; i < Tr.NbrL; i++){
            fou <<  i+1 << endl;
            fou <<  Tr.B[i+2].dx << endl;
            nd = Tr.B[i+2].pts;
            fou <<  nd << endl;
            fou <<  "Qave" << endl;
            fou <<  Tr.B[i+2].Qave  << endl;
            fou <<  "FFR" << endl;

             for (long j = 0; j < Tr.B[i+2].pts; j++){

                Pd = Tr.B[i+2].Pave[j];
                FFR = Pd/Pa;
                fou <<  FFR << endl;

             }
        }

        fou <<  "RCA" << endl;
        for (long i = 0; i < Tr.NbrR; i++){
            fou <<  i+1 << endl;
            fou <<  Tr.B[i + 2 + Tr.NbrL].dx << endl;
            nd = Tr.B[i + 2 + Tr.NbrL].pts;
            fou <<  nd << endl;
            fou <<  "Qave" << endl;
            fou <<  Tr.B[i + 2 + Tr.NbrL].Qave  << endl;
            fou <<  "FFR" << endl;

             for (long j = 0; j < Tr.B[i + 2 + Tr.NbrL].pts; j++){

                Pd = Tr.B[i + 2 + Tr.NbrL].Pave[j];
                FFR = Pd/Pa;
                fou <<  FFR << endl;

             }
        }

        fou.close();

	}
};
