//########################################################################
//  ������ ����������� ������-��������� ����������
//  ������-���������� ���������� ���������, ������� ��� ����������
//  ������� ������������ ����� ������
//########################################################################
#include <string>
#include <vector>
#include "CalcUtils.h"

namespace TreeDataDepend {
	using namespace std;
	using namespace CalcUtils;

	const static long IT_BLOOD = 1;
	const static long IT_SUBSTANCE = 2;
	const static long GT_LUMPED = 1;
	const static long GT_DISTRIBUTED = 2;

	//________________________________________________________________________
	//________________________________________________________________________
	//     ������-��������� ��������� �����
	//________________________________________________________________________
	class TDVetvParameter {
	public:
		double R , Rinit;   // ������������� (���������������� � ������ ����� ������)
		double delta;    // ����������� ����������� � ������
		double c , cnewQ , cnewP;      // �������� ��������� ����� � �������� �����������
		double Ps;      // ����. ��� ������� ������� ���������
		double Smin , Smax;  // ������������ � ����������� �������� ����������� ������� ������
		double S0 , S0new; // ���������� ������� � �����
		double Par , Pgr;   // ������������ (Par) � ����������� (Pgr) ��������, � ������� ����� ������ ���� 0
		double d_start;
		long gravity , pmaxflag; // ����������
		double Qave_curr , Rave_curr , Pave_curr , RadiusAve_curr;
		double Qave_prev , Rave_prev , Pave_prev , RadiusAve_prev;
		double Qave0 , Qave_next , Pave0 , Pave_next;
		double Save_curr , Save_next , Save_prev , Save0;
		double *Pmax , *Pmaxcurr;
		double MDP , PD_aor_0 , PD_aor_EECP , MSP , PS_aor_0 , PS_aor_EECP , EE , V_0 , V_comp; //// EECP
		long EEflag;
		vector< vector<double> > conc , conc_old; // ������������ �������

		vector<string> fnameC; // ������������ ������� ������������ � ��� �����
		vector<string> fnameVar; //string fnameS , fnameU , fnameP;
		vector<string> fnameTD; //*fnamePmax , *fnamePave , *fnameQave , *fnameMDP , *fnameMSP , *fnameEE;
		string strID;
	};

	//________________________________________________________________________
	//     ���������� ������������� ������
	//________________________________________________________________________
	class TDAlveolar {
	public:
		double A , I , Ig , R , Rg , C , Y [ 4 ];
	};

	//________________________________________________________________________
	//     ������-��������� ��������� ����
	//________________________________________________________________________
	class TDUzelParameter {
	public:
		long TextID;         // ����� ������ ����������� ������� ����
		long Kext_sz;         // ���������� ������� �����
		vector<long> KextID; // ������ ������� ����� (�� ������ �������� ������)
		TDAlveolar *alv;
		long NgasExch;
		vector<long> gasExchKnt;
		// LYMPHATIC
		double LymJVdst;
	};

	//________________________________________________________________________
	//     ������� ����������� �� �����
	//________________________________________________________________________
	class TDExternalImpact {
	public:
		long BrID;     // ����� ������
		long IType;     // ��� ����������� (���� / ������ ����� / ������������ �������)
		long GType;     // ��� ��������� �����������: 1 - ��������, 2 - ������������ �� ���� ����� ������
		double len;       // ���� ��� ��������� ��������, �� ����� ����������� (� % ����� ������, �� ��� ������)
		long SubstID;    // ����� �������� (���� ���� / ������ ������������)
		double alfa;       // ����������� ������������ ����������� (�������� ����� / ������� � �.�.)
		double Tstart , Tend;   // ������� ������ � ��������� ����������� (����������)
	};
}
