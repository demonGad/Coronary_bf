#include "TDGranuslov.h"

namespace Granuslov {

    using namespace CalcUtils;
    using namespace TaskData;
    using namespace TDGranuslov;


    /*
    Calculate function F in a X for Newton's method, Bernoulli
    */

    vector<double> FUZ_Ber (Zadacha& Z, vector<double> Xs, vector<double> Xu, vector<Vetv*> brin, vector<Vetv*> brou, long N){

        long szin, szou, j;

        vector<double> alf(N), bet(N), e(N), F(N), P(N), out(3);

        double rho;

        rho = 1.0;

        szin = brin.size();
        szou = brou.size();

        for (long i = 0; i < szin; i++){

            e[i] = 1;
            F[i] = 0;
            tie(alf[i],bet[i]) = IncomingCompatibilityCoeffs(Z, *(brin[i]));

            out = (*(brin[i])).URSOB(Xs[i],Xu[i]);
            P[i] = out[0];
        }

        for (long i = 0; i < szou; i++){

            j = i + szin;
            e[j] = - 1;
            F[j] = 0;
            tie(alf[j],bet[j]) = OutgoingCompatibilityCoeffs(Z, *(brou[i]));

            out = (*(brou[i])).URSOB(Xs[j],Xu[j]);
            P[j] = out[0];
        }

        F[0] = e[0]*(alf[0]*Xs[0] + bet[0])*Xs[0];

        for (long i = 1; i < N; i++){

            F[i] = P[i-1]/rho - P[i]/rho + 0.5*pow(alf[i-1]*Xs[i-1] + bet[i-1],2) - 0.5*pow(alf[i]*Xs[i] + bet[i],2);
            F[0] = F[0] + e[i]*(alf[i]*Xs[i] + bet[i])*Xs[i];
        }

        return F;
    }


    /*
    Calculate Jacobian J in a X for Newton's method, Bernoulli
    */

    Matrix TDYacobian_Ber(Zadacha& Z, vector<double> Xs, vector<double> Xu, vector<Vetv*> brin, vector<Vetv*> brou, long N){

        long szin, szou, j;

        vector<double> alf(N), bet(N), e(N), F(N), P(N), Ps(N), out(3);

        Matrix YAC(N,N);

        double rho;

        rho = 1.0;

        szin = brin.size();
        szou = brou.size();

        for (long i = 0; i < szin; i++){

            e[i] = 1;
            F[i] = 0;
            tie(alf[i],bet[i]) = IncomingCompatibilityCoeffs(Z, *(brin[i]));

            out = (*(brin[i])).URSOB(Xs[i],Xu[i]);
            P[i] = out[0];
            Ps[i] = out[1];
        }

        for (long i = 0; i < szou; i++){

            j = i + szin;
            e[j] = - 1;
            F[j] = 0;
            tie(alf[j],bet[j]) = OutgoingCompatibilityCoeffs(Z, *(brou[i]));

            out = (*(brou[i])).URSOB(Xs[j],Xu[j]);
            P[j] = out[0];
            Ps[j] = out[1];
        }

        YAC.Zero();

        for (long i = 0; i < N; i++)
            YAC(0,i) = e[i]*(2*alf[i]*Xs[i]+ bet[i]);

        for (long i = 1; i < N; i++){

            YAC(i,i-1) = Ps[i-1]/rho + (alf[i-1]*Xs[i-1] + bet[i-1])*alf[i-1];
            YAC(i,i) = (-1)*(Ps[i]/rho + (alf[i]*Xs[i] + bet[i])*alf[i]);
        }

        return YAC;
    }


    double vecNorm (vector<double> vec){

        double cur_norm;

        cur_norm = 0;

        for (long i = 0; i < vec.size(); i++)
            if (fabs(vec[i]) > cur_norm)
                cur_norm = abs(vec[i]);

        return cur_norm;

    }

    vector<double> vecSum (vector<double> x1, vector<double> x2, int N){

        vector<double> res(N,0.0);

        for (long i = 0; i < N; i++)
            res[i] = x1[i] + x2[i];

        return res;
    }

    vector<double> vecSc (vector<double> x1, double a, int N){

        vector<double> res(N,0.0);

        for (long i = 0; i < N; i++)
            res[i] = x1[i]*a;

        return res;
    }

    /*
    Calculate velocities at the end of Newton's method iteration
    */
    vector<double> TDNewtonPostproc (Zadacha& Z, vector<double> Xs, vector<Vetv*> brin, vector<Vetv*> brou, long N){

        vector<double> U(N);
        double alfa, beta;
        long szin, szou, j;

        szin = brin.size();
        szou = brou.size();

        for (long i = 0; i < szin; i++){

            tie(alfa,beta) = IncomingCompatibilityCoeffs(Z, *(brin[i]));

            U[i] = alfa*Xs[i] + beta;
        }

        for (long i = 0; i < szou; i++){

            j = i + szin;
            tie(alfa,beta) = OutgoingCompatibilityCoeffs(Z, *(brou[i]));

            U[j] = alfa*Xs[j] + beta;
        }

        return U;

    }


    void CalculateCommonKnot(Zadacha& Z, long kID, vector<Vetv*> brin, vector<Vetv*> brou, long N){

        long szin, szou, Nmax, i, p;

        double Det, X0_norm, F_norm, F2_norm, F_2_norm, eps;

        Matrix YAC(N,N), YAC_inv(N,N);

        vector<double> F(N), X0(N), Xs(N), Xu(N), X0_cur(N), F2(N), F_2(N);

        szin = brin.size();
        szou = brou.size();

        Nmax = 1000;

        for ( long i = 0; i < szin; i++ ){

           // Xs[i] = (*(brin[i])).VBO[0][(*(brin[i])).pts - 2];
            //Xu[i] = (*(brin[i])).VBO[1][(*(brin[i])).pts - 2];

            Xs[i] = (*(brin[i])).VBO[0][(*(brin[i])).pts - 1];
            Xu[i] = (*(brin[i])).VBO[1][(*(brin[i])).pts - 1];
        }

        for ( long i = 0; i < szou; i++ ){

            //Xs[szin + i] = (*(brou[i])).VBO[0][1];
            //Xu[szin + i] = (*(brou[i])).VBO[1][1];

            Xs[szin + i] = (*(brou[i])).VBO[0][0];
            Xu[szin + i] = (*(brou[i])).VBO[1][0];
        }

        i = 1;
        X0_norm = 1.0;
        eps = 0.00001;


        //cout << "Here before while  " <<  T <<  "  id = "<<  kID <<endl;
        while ((i < Nmax)&&(X0_norm > eps)){

            F = FUZ_Ber(Z, Xs, Xu, brin, brou, N);


            YAC = TDYacobian_Ber(Z, Xs, Xu, brin, brou, N);

            YAC_inv = YAC.InvMatrix();

            for (long j = 0; j < N; j++)
                F[j] = (-1)*F[j];

            X0 = YAC_inv.MulVr(F);

            p = 0;
            //cout << "Here before p  " <<  T <<  "  id = "<<  kID <<endl;
            while (p < 10){

                F = FUZ_Ber(Z, vecSum(Xs,X0,N), Xu, brin, brou, N);

                F2 = FUZ_Ber(Z, vecSum(Xs,vecSc(X0,2,N),N), Xu, brin, brou, N);

                F_2 = FUZ_Ber(Z, vecSum(Xs,vecSc(X0,0.5,N),N), Xu, brin, brou, N);

                F_norm = vecNorm(F);
                F2_norm = vecNorm(F2);
                F_2_norm = vecNorm(F_2);

                if ((F2_norm <= F_norm)&&(F2_norm <= F_2_norm)){
                    X0 = vecSc(X0,2,N);
                    p++;
                }
                if ((F_2_norm < F_norm)&&(F_2_norm < F2_norm)){
                    X0 = vecSc(X0,0.5,N);
                    p++;
                }
                if ((F_norm <= F_2_norm)&&(F_norm <= F2_norm)){
                    p = 10;
                }

            }
            //cout << "Here after p  " <<  T <<  "  id = "<<  kID <<endl;


            for (long j = 0; j < N; j++)
                Xs[j] = Xs[j] + X0[j];

            X0_norm = 0;

            for (long j = 0; j < N; j++)
                X0_norm = X0_norm + X0[j]*X0[j];

            X0_norm = sqrt(X0_norm);

            Xu = TDNewtonPostproc(Z,Xs, brin, brou, N);

            i++;

        }

        //cout << "Here after while  " <<  T <<  "  id = "<<  kID << "  n_iter = "<<  i <<endl;

        if (i > Nmax)
            cout << "Warning [CalculateCommonKnot]: Iteration Limit ; kID " << kID << endl;

        /*if (T > 0.001){
        //if (kID == 2){
            for ( long i = 0; i < N; i++ ){

                cout << i << "   Xs = " << Xs[i] <<"   Xu = " << Xu[i] << endl;
            }
            cout << "  T =  "  << T << endl;
        }*/


        // save results
        for ( long i = 0; i < szin; i++ ){

            (*(brin[i])).VB[0][(*(brin[i])).pts - 1] = Xs[i];
            (*(brin[i])).VB[1][(*(brin[i])).pts - 1] = Xu[i];

            if (Xs[i] < 0){
                cout << "Xs[i] < 0   BrId = "  <<  (*(brin[i])).ID << "; kID =  " << kID << endl;
            }

            if ((Xs[i]/ (*(brin[i])).TD.S0 ) < 0.3) {
                cout << "Xs[i] is small   BrId = "  <<  (*(brin[i])).ID << "; kID =  " << kID << endl;
            }
        }

        for ( long i = 0; i < szou; i++ ){

            (*(brou[i])).VB[0][0] = Xs[szin + i];
            (*(brou[i])).VB[1][0] = Xu[szin + i];
        }


    }

    inline void Grtoch(Zadacha& Z , Derevo& Tr , Uzel& kn){


        if ((kn.IG == FLOW)&&((kn.Nin + kn.Nou) == 1)){
            //cout << "Here Bound  " <<  T <<  "  id = "<<  kn.ID <<endl;
            TDGrtoch(Z, Tr, kn); // inner or outer knot
        }
        else if ((kn.Nin + kn.Nou) > 1)
        {

            //CallID is removed
            //cout << "Here Inner  " <<  T <<  "  id = "<<  kn.ID <<endl;
            CalculateCommonKnot(Z, kn.ID, kn.Bin, kn.Bou, kn.Nou + kn.Nin); // simplified, kn.ID - for error messages

            /* Debug

            cout << "kn.ID  " << kn.ID << endl;

            for (long i = 0; i < kn.Bin.size(); i++ ){

                    cout << "IDin  " << (*kn.Bin[i]).ID << endl;
                    cout << "S  " << (*(kn.Bin[i])).VB[0][(*(kn.Bin[i])).pts - 1];
                    cout << "   U  " << (*(kn.Bin[i])).VB[1][(*(kn.Bin[i])).pts - 1] << endl;
            }

            for (long i = 0; i < kn.Bou.size(); i++ ){

                    cout << "IDou  " << (*kn.Bou[i]).ID << endl;
                    cout << "S  " << (*(kn.Bou[i])).VB[0][0];
                    cout << "   U  " << (*(kn.Bou[i])).VB[1][0] << endl;
            }

            cout << "kn.ID  " << kn.ID << endl;
            */


        }
    }
};
