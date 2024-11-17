

namespace TDGranuslov {

    using namespace CalcUtils;
    using namespace TaskData;
    using namespace Subroutines;

    double getInFlow (Zadacha& Z){

            double Qin, Tc, Tcur, Ts, Tf;

            if (T > (Z.T_last_Hbeat + Z.h_period_curr)){

                // new heart cycle

                Z.T_last_Hbeat = Z.T_last_Hbeat + Z.h_period_curr;
                Z.N_heart_cycles += 1;
            }

            Tc = 60/Z.HR;
            Tcur = (T - Z.T_last_Hbeat)/Tc;
            Ts = 0;
            Tf = 0.3;

            if ((Tcur > Ts)&&(Tcur < Tf))
                Qin = (Z.SV*PI/(2*(Tf - Ts)))*sin(PI*(Tcur - Ts)/(Tf - Ts));
            else
                Qin = 0;

            return Qin;
    }

    vector<double> TDFlowToSU (Zadacha& Z , Derevo& Tr , Uzel& kn, double Qin){

        vector<double> V(Z.Cor);
        double alfa, beta;

        if (kn.Nou == 1){ // beginning of a branch

            tie(alfa, beta) = OutgoingCompatibilityCoeffs(Z,*(kn.Bou[0]));

            V[1] = (beta + sqrt(beta*beta + 4*alfa*Qin))/2;
            V[0] = (- beta + sqrt(beta*beta + 4*alfa*Qin))/(2*alfa);
        }

        if (kn.Nin == 1){

            tie(alfa, beta) = IncomingCompatibilityCoeffs(Z,*(kn.Bin[0]));

            V[1] = (beta - sqrt(beta*beta + 4*alfa*Qin))/2;
            V[0] = (- beta - sqrt(beta*beta + 4*alfa*Qin))/(2*alfa);

        }

        return V;
    }



    // Calculate function and derivative for Newton's method

    tuple<double, double> FuncNewtS(double Scur, Vetv& br, double alfa, double beta, double Pout){

        double Ps, Pcur, Func, FuncS;
        vector<double> out(3);

        out = br.URSOB(Scur, br.VBO[1][br.pts - 1]);
        Pcur = out[0];
        Ps = out[1];

        if (Pcur < 0)
            cout << "Warning [FuncNewtS]: P < 0 ; brID " << br.ID << endl;

        Func = alfa*Scur*Scur + beta*Scur - (Pcur - Pout)/br.TD.R;
        FuncS = 2*alfa*Scur + beta - Ps/br.TD.R;

        return {Func, FuncS};

    }

    // Newton's method for outlet Q = (P - Pveins)/R
    vector<double> CalcRtoSUnewt (Zadacha& Z , Derevo& Tr , Vetv& br){

        vector<double> V(Z.Cor);
        double Pout, S_norm, Scur, Snew, alfa, beta, Func, FuncS, Pcor, eps;
        long i, Nmax;

        i = 0;
        Nmax = 10000;
        eps = 0.001;

        //Pout = Z.Pveins*1333.2;

        Pout = Pout = 10665;

        Pcor = 0; // external Pressure, we assume it to be 0

        tie(alfa,beta) = IncomingCompatibilityCoeffs(Z,br);

        Scur = ( -beta - sqrt(beta * beta ) ) / (2. * alfa);
        Snew = br.TD.S0*(1 + log(1 + (Pout - Pcor)/(br.TD.c*br.TD.c)));

        /*cout << "BrID: " << br.ID << endl;
        cout << "Scur: " << Scur << "   Snew:  "<< Snew <<  endl;
        cout << "alfa: " << alfa << "   beta:  "<< beta <<  endl;
        cout << endl;*/


        if ((Scur <= 0)||(Snew <=0))
            cout << "Warning [CalcRtoSUnewt]: S <= 0 ; brID " << br.ID << endl;

        if (Scur < Snew)
            Scur = Snew;

        S_norm = 1;

        while ((i < Nmax)&&(S_norm > eps)){

            tie(Func,FuncS) = FuncNewtS(Scur,br,alfa,beta,Pout);

            Snew = Scur - (Func/FuncS);

            if (abs(Scur) > 0.00000001)
                S_norm = abs((Snew - Scur)/Scur);
            else{
                cout << "Warning [CalcRtoSUnewt]: Scur = 0 ; brID " << br.ID << endl;
                S_norm = 0;
            }

            i++;
            Scur = Snew;
        }

        if (i > Nmax)
            cout << "Warning [CalcRtoSUnewt]: Iteration Limit ; brID " << br.ID << endl;

        V[0] = Snew;
        V[1] = alfa*Snew + beta;

        return V;
    }

    //calculate boundary knots
    void TDGrtoch(Zadacha& Z , Derevo& Tr , Uzel& kn){

        long idx1,idx2;
        double Qin, alf, bet, Tc, Qin_coeff, Pout, Rcor, Tcur, Ts, Tf;
        vector<double> V(Z.Cor);

        Vetv *brp = &(Tr.B[0]);

        // arterial network

        if ((Tr.ID == PULMART)||(Tr.ID == SYSART)){

            if ((kn.IG == FLOW)&&(kn.Nou == 1)){

                //aorta
                brp = kn.Bou[0];

                idx1 = 0;
                idx2 = 1;

                Tc = 60/Z.HR;
                //if (Tc != 1.0)
                //    cout << "Tc is odd, Tc = " << Tc << endl;
                Qin = getInFlow(Z);

                //if (T > 5.0)
                //    Qin = 0;

                V = TDFlowToSU(Z, Tr, kn, Qin);

                if (V[0] < 0){
                    cout << "V[0] < 0 !!! FlowToSU " <<endl;

                }
            }
            if ((kn.IG == FLOW)&&(kn.Nin == 1)){

                //terminal artery

                brp = kn.Bin[0];

                idx1 = (*brp).pts - 1;
                idx2 = (*brp).pts - 2;

                Tc = 60/Z.HR;
                if (Tc != 1.0)
                    cout << "Tc is odd, Tc = " << Tc << endl;
                Tcur = (T - Z.T_last_Hbeat)/Tc;
                Ts = 0;
                Tf = 0.3;

                Rcor = 0;
                if ((Tcur > Ts)&&(Tcur < Tf))
                    Rcor = sin(PI*(Tcur - Ts)/(Tf - Ts));

                if ((*brp).group == 1) //RCA
                    (*brp).TD.R = (*brp).TD.Rinit*(1 + 5*Rcor);

                if ((*brp).group == 2) //LCA
                    (*brp).TD.R = (*brp).TD.Rinit*(1 + 10*Rcor);

                V = CalcRtoSUnewt(Z, Tr, *brp);
            }

            (*brp).VB[0][idx1] = V[0];
            (*brp).VB[1][idx1] = V[1];

           /* cout << "BrID: " << (*brp).ID << endl;
            cout << "Brlen: " << (*brp).len << "   BrD:  "<< (*brp).width <<  endl;
            cout << "V: " << V[0] << "  " << V[1] << endl;
            cout << endl;*/

        }

    }
};
