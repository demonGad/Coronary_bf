#include <string>
#include <cmath>  

namespace TDGlobals {
	using namespace std;

	string root_out , GPLOT_root_out;
	long isFirstTime;
	double HeartPeriod;
	long NBeats_prev;
	double T_last_Hbeat; //����� ���������� ���������� ���������� �����
	double h_period_curr; //����� �������� ���������� �����
	long N_heart_cycles; //����� ��������� ������ (�����������) ; ���� ����� ��������� ��� ���������� � ���� ������ "������" 
	double rho_air;
	double timeForTrack;

	// LYMPHATIC SYSTEM
	// Pressure equation P = teta(S-S0) + Ampl*sin(omega(X - vel*T))
	double teta , Pin , Pout , Ampl , vel , omega;
	inline double LymPofS (double S, double S0, double x, double t) {
		return teta * ( S - S0 ) + Ampl * sin( omega * ( x - vel * t ) );
	}
	inline double LymdPofdS () {
		return teta;
	}
	inline double LymSofP (double P, double S0, double x, double t ) {
		return S0 + 1.0 / teta * ( P - Ampl * sin( omega * ( x - vel * t ) ) );
	}
}