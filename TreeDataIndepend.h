//########################################################################
//	������ ����������� ������-����������� ��������� ������ ��� ���� (�����)
// Uzel   - ���� (����� ���������, ��������� �����)
// Vetv   - ����� (�����)
// Derevo - �� ���� ���� (����), ������ - ������� ������, �������� ������
// ������, �����, ������� ����������� � �.�., � ����� ����� ������ �
// ����������� ������ � �����
//########################################################################
#include "TreeDataDepend.h"
#include <vector>
//#include <math.h>
namespace TreeDataIndepend {

	using namespace TreeDataDepend;
	//________________________________________________________________________
	//					����
	//	�		- ����� �� ������������ ���������
	//	ID	- ����� ����
	// Nou	- ���-�� ��������� ������
	// Nin	- ���-�� �������� ������
	// IG	- ��� ���� (���������, ����, ����� � �.�.)
	// TD	- ������-��������� ����� ���������� ����
	// ���������� ��� ��������� ����� � ���� (�������� ��� ���������)
	//________________________________________________________________________
	class Vetv;
	class Uzel {
	public:
		Dot C;
		long ID , Nou , Nin , IG;
		vector<Vetv*> Bou;
		vector<Vetv*> Bin;
		//bool LYMPHNODE;
		TDUzelParameter TD;

	};

	//________________________________________________________________________
	//					�����
	//	ID				- ����� �����
	// pts				- ���������� ����� �� �����
	// dx				- ��� �� �����
	// len				- ����� �����
	//	width			- ������ �����
	//	Kn1, Kn2	- ��������� �� ����������� � ����� ����; Kn1 - ������, Kn2 - �����
	//	VB				- �������� ������� ���������� �� �����
	// VBO				- �������� ������� ���������� �� ����� �� ���������� ����
	// TD				- ������-��������� ����� ���������� �����
	// gravity           - �������� �� ���������� (1 ��� 0)
	//________________________________________________________________________
	class Vetv {
	public:
		long ID , pts , myTreeID , InvertPoints, group, stenType;
		Uzel *Kn1 , *Kn2;
		vector< vector<double> > VB, VBO;
		vector<double>  Pave;
		double len , width , dx, Qave;
		TDVetvParameter TD;


		// S - cross-section, u - velocity
		vector<double> URSOB(double S,double u) const {

		    vector<double> output(3);       // [0] - pressure, [1] - derivative Pr, [2] - zvuk
		    double rho_w,rho;                   //density of a wall

		    rho_w = 1.0;
		    rho = 1.0;

		    if (S > TD.S0){

                output[0] = TD.c*TD.c*rho_w*(exp(S/TD.S0 - 1) - 1);
                output[1] = TD.c*TD.c*rho_w*exp(S/TD.S0 - 1)/TD.S0;

		    } else {

                output[0] = TD.c*TD.c*rho_w*log(abs(S/TD.S0));
                output[1] = TD.c*TD.c*rho_w/S;
		    }

		    output[2] = sqrt(S*output[1]/rho);

		    if ( output[2] < 0.00001){
                output[2] = 0.00001;
                cout << "Warning! PWV is too small" << endl;
		    }

			return output;
		}


		vector<double> FPRCH(double S,double u, long i) const {

		    vector<double> FPR(2);

		    FPR[0] = 0;
		    FPR[1] = (-8.0)*(3.1415926)*(0.04)*u/S;

		    return FPR;
		}
	};

	void printVetv ( Vetv& B ) {
		cout << "ID = " << B.ID << " pts = " << B.pts << " myTreeID = " << B.myTreeID << endl;
		cout << "len = " << B.len << " width = " << B.width << " dx = " << B.dx << endl;
		cout << "Kn1 = " << B.Kn1->ID << " Kn2 = " << B.Kn2->ID << endl;
		cout << "--------------------------------------" << endl;
	}
	void printUzel ( Uzel& K ) {
			cout << "ID = " << K.ID << " IG = " << K.IG <<  endl;
			cout << "Nin = " << K.Nin << " Nou = " << K.Nou << endl;
			cout << "Bin : ";
			for (int i = 0; i < K.Bin.size(); ++i) cout << K.Bin[i]->ID << " ";
			cout << endl <<  "Bou : ";
			for (int i = 0; i < K.Bou.size(); ++i) cout << K.Bou[i]->ID << " ";
			cout << endl << "--------------------------------------" << endl;
		}

	//________________________________________________________________________
	//					������
	//	ID	- ����� ������
	// Nbr	- ���������� ������
	// Nkn	- ���������� �����
	//	B	- ������ ������ ������
	//	K	- ������ ����� ������
	//________________________________________________________________________
	class Derevo {
	public:
		long ID , Nbr , Nkn , Nimp , Norg, NbrL, NbrR, NknL, NknR, Nlad, Nlcx, Nrca, iLad, iLcx, iRca, iCa;
		string treefilename , dirname;
		string knotfilename , TDknotfilename , TDknotExternalFilename;
		string branchfilename , TDbranchfilename;
		string TDImpactfilename;
		string inputdata;
		string FFRfile, FFRfullfile;
		Vetv *B;
		Uzel *K;
		vector<TDExternalImpact> I;
	};

	// ������ ����� ������ ������ ��� ����������
	class MultiKnotTreeGroup {
	public:
		long TrID;
		long sz;
		vector<Uzel> KntLst;
	};

	// ����������  - ������������ ������� ���������� �������� ����� ������ ������ (��� ������ � ���� �� �����) � ����
	class MultiKnot {
	public:
		long sz;
		long szin;
		long szou;
		vector<MultiKnotTreeGroup> GrpLst;
	};

	// ������ �����������
	class MultiKnotList {
	public:
		long sz;
		vector<MultiKnot> Lst;
		string multiKnotsfilename;
	};
};
