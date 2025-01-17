#include "TreeInitTD.h"

namespace TreeInit {
	using namespace TreeInitTD;

	inline void TreeInitialization( Derevo& Tr , Zadacha& Z ) {
		Uzel *knt1 , *knt2;
		string NID;
		long tmp_i, addBr, addKn;
		double tmp_d;

		long ID , KnID1 , KnID2;


        Tr.dirname = "tree1" + trim( Adjustl( NID ) );
		/* not used in coronary
		{
			char wrotmnenogi [ 10 ];
			sprintf( wrotmnenogi , "%ld" , Tr.ID );
			NID = string( wrotmnenogi );
		}
		if ( Tr.ID == BRONCHIAL ) {
			NID = "bronch";
			//InitAlveolar();
		}
		// trim == long to string //
        Tr.dirname = "tree" + trim( Adjustl( NID ) );
		Tr.treefilename = trim( SharedDirectory ) + "tree" + trim( Adjustl( NID ) ) + slash + "tree.tre";
		Tr.knotfilename = trim( SharedDirectory ) + "tree" + trim( Adjustl( NID ) ) + slash + "knot.tre";
		
		Tr.TDknotExternalFilename = trim( SharedDirectory ) + "tree" + trim( Adjustl( NID ) ) + slash + "tdknotexternal.tre";
		Tr.TDbranchfilename = trim( SharedDirectory ) + "tree" + trim( Adjustl( NID ) ) + slash + "tdbranch.tre";
		Tr.TDImpactfilename = trim( SharedDirectory ) + "tree" + trim( Adjustl( NID ) ) + slash + "tdimpact.tre";*/
		Tr.inputdata =  trim( SharedDirectory ) + "input_data.tre";
		Tr.FFRfile =  trim( SharedDirectory ) + "FFR.tre";
		Tr.FFRfullfile =  trim( SharedDirectory ) + "FFRfull.tre";
		Tr.branchfilename = trim( SharedDirectory ) + "tree1" + trim( Adjustl( NID ) ) + slash + "branch.tre";

        ofstream fou;
        fou.open(Tr.FFRfile);
        fou << (-1) << endl;
        fou.close();

        fou.open(Tr.FFRfullfile);
        fou << (-1) << endl;
        fou.close();


		cout << "Intialization of tree " << endl;
		ifstream fin( Tr.branchfilename, ifstream::in );
		fin >> tmp_i;
		Tr.Nkn = tmp_i + 1;
		fin >> tmp_i;
		Tr.Nbr = tmp_i + 2;

		Tr.B = new  Vetv [ Tr.Nbr ];
		Tr.K = new  Uzel [ Tr.Nkn ];

		//----------------------------------- Reading Knots ----------------------------------

		//LCA and Aorta
		fin >> Tr.NknL;

        addKn = 0;																							//UPDATE
		for ( long i = 1; i < Tr.NknL; i++ ) {
			fin >> ID;
			ID = ID + addKn;
			Tr.K [ ID - 1 ].ID = ID;
			fin >> Tr.K [ ID - 1 ].C.X;
			fin >> Tr.K [ ID - 1 ].C.Y;	//
			fin >> Tr.K [ ID - 1 ].C.Z;	//
			fin >> Tr.K [ ID - 1 ].IG; //
		}

		//RCA
		fin >> Tr.NknR;

		fin >> tmp_i;  // skip first node, we have it
		fin >> tmp_d;
		fin >> tmp_d;
		fin >> tmp_d;
		fin >> tmp_i;

		addKn = Tr.NknL + 1;

		for ( long i = 1; i < Tr.NknR; i++ ) {
			fin >> ID;
			ID = ID + addKn;
			Tr.K [ ID - 1 ].ID = ID;
			fin >> Tr.K [ ID - 1 ].C.X;
			fin >> Tr.K [ ID - 1 ].C.Y;	//
			fin >> Tr.K [ ID - 1 ].C.Z;	//
			fin >> Tr.K [ ID - 1 ].IG; //
		}


		//----------------------------------- Reading Branches  --------------------------------


		//Aorta																										//UPDATE
        //LCA
		fin >> Tr.NbrL;

        addBr = 0;
        addKn = 0;
		for ( long i = 0; i < Tr.NbrL; i++ ) {
			fin >> ID;
			ID = ID + addBr;
			Tr.B [ ID - 1 ].ID = ID;
			Tr.B [ ID - 1 ].myTreeID = Tr.ID;
			Tr.B [ ID - 1 ].InvertPoints = 0;

			fin >> KnID1 >> KnID2;		//
			fin >> Tr.B [ ID - 1 ].len;		//
			fin >> Tr.B [ ID - 1 ].width;	//
			fin >> Tr.B [ ID - 1 ].pts;		//


			if ( ( Tr.B [ ID - 1].len > 15.0 ) && ( Tr.B [ ID - 1 ].pts < 15 ) ) {
				Tr.B [ ID - 1 ].pts = 15;
			}

			if ( Tr.B [ ID - 1].len < 1.2 )
                Tr.B [ ID - 1 ].pts = 8;

            if ( Tr.B [ ID - 1].len < 0.7 )
                Tr.B [ ID - 1 ].pts = 5;

            if ( Tr.B [ ID - 1].len < 0.35 )
                Tr.B [ ID - 1 ].pts = 3;


			Tr.B [ ID - 1 ].Kn1 = &Tr.K [ KnID1 - 1 + addKn ];
			Tr.B [ ID - 1 ].Kn2 = &Tr.K [ KnID2 - 1 + addKn ];

			if ( Tr.B [ ID - 1 ].pts > 1 ) {
				Tr.B [ ID - 1 ].dx = Tr.B [ ID - 1 ].len / double( Tr.B [ ID - 1 ].pts - 1 );
			}
			else {
				cout << "Error! Wrong points number at (tree, branch): " << Tr.ID << Tr.B [ ID - 1 ].ID << endl;
				Tr.B [ ID - 1 ].dx = Tr.B [ ID - 1 ].len;
			}

			if (Tr.K[KnID2 - 1 + addKn].IG == 2){
                Tr.B[ID - 1].group = 2;  // LCA terminal
			}else{
                Tr.B[ID - 1].group = 4;  // LCA non-terminal
			}

			 Tr.B[ID - 1].stenType = 0;
			 Tr.B[ID - 1].Pave.resize( Tr.B [ ID - 1 ].pts);                                             ///Поломка
			 Tr.B[ID - 1].Qave = 0;
			 for ( long j = 0; j < Tr.B [ ID - 1 ].pts; j++ ){
                Tr.B[ID - 1].Pave[j] = 0;
            }


			{
				Tr.B [ ID - 1 ].VB.resize( Z.Cor );
				Tr.B [ ID - 1 ].VBO.resize( Z.Cor );
				for ( long j = 0; j < Z.Cor; j++ ) {
					Tr.B [ ID - 1 ].VB [ j ].resize( Tr.B [ ID - 1 ].pts );
					Tr.B [ ID - 1 ].VBO [ j ].resize( Tr.B [ ID - 1 ].pts );
				}
			}
		} // LCA for-loop


		Tr.B[0].Kn1 = &Tr.K[0]; //connect LCA root to aorta										 //UPDATE

		//RCA

		fin >> Tr.NbrR;



		addBr = Tr.NbrL;																		//UPDATE
        addKn = Tr.NknL + 1;

		for ( long i = 0; i < Tr.NbrR; i++ ) {
			fin >> ID;
			ID = ID + addBr;
			Tr.B [ ID - 1 ].ID = ID;
			Tr.B [ ID - 1 ].myTreeID = Tr.ID;
			Tr.B [ ID - 1 ].InvertPoints = 0;

			fin >> KnID1 >> KnID2;		//
			fin >> Tr.B [ ID - 1 ].len;		//
			fin >> Tr.B [ ID - 1 ].width;	//
			fin >> Tr.B [ ID - 1 ].pts;		//

			if ( ( Tr.B [ ID - 1].len > 15.0 ) && ( Tr.B [ ID - 1 ].pts < 15 ) ) {
				Tr.B [ ID - 1 ].pts = 15;
			}

			if ( Tr.B [ ID - 1].len < 1.2 )
                Tr.B [ ID - 1 ].pts = 8;

            if ( Tr.B [ ID - 1].len < 0.7 )
                Tr.B [ ID - 1 ].pts = 5;

            if ( Tr.B [ ID - 1].len < 0.35 )
                Tr.B [ ID - 1 ].pts = 3;


			Tr.B [ ID - 1 ].Kn1 = &Tr.K [ KnID1 - 1 + addKn ];
			Tr.B [ ID - 1 ].Kn2 = &Tr.K [ KnID2 - 1 + addKn ];

			if ( Tr.B [ ID - 1 ].pts > 1 ) {
				Tr.B [ ID - 1 ].dx = Tr.B [ ID - 1 ].len / double( Tr.B [ ID - 1 ].pts - 1 );
			}
			else {
				cout << "Error! Wrong points number at (tree, branch): " << Tr.ID << Tr.B [ ID - 1 ].ID << endl;
				Tr.B [ ID - 1 ].dx = Tr.B [ ID - 1 ].len;
			}

			if (Tr.K[KnID2 - 1 + addKn].IG == 2){
                Tr.B[ID - 1].group = 1;  // RCA terminal
			}else{
                Tr.B[ID - 1].group = 3;  // RCA non-terminal
			}

			 Tr.B[ID - 1].stenType = 0;

			 Tr.B[ID - 1].Pave.resize( Tr.B [ ID - 1 ].pts);
			 Tr.B[ID - 1].Qave = 0;
			 for ( long j = 0; j < Tr.B [ ID - 1 ].pts; j++ ){
                Tr.B[ID - 1].Pave[j] = 0;
			 }

			{
				Tr.B [ ID - 1 ].VB.resize( Z.Cor );
				Tr.B [ ID - 1 ].VBO.resize( Z.Cor );
				for ( long j = 0; j < Z.Cor; j++ ) {
					Tr.B [ ID - 1 ].VB [ j ].resize( Tr.B [ ID - 1 ].pts );
					Tr.B [ ID - 1 ].VBO [ j ].resize( Tr.B [ ID - 1 ].pts );
				}
			}
		} // RCA for-loop

        Tr.B[Tr.NbrL + 1 - 1].Kn1 = &Tr.K[0]; // connect to aorta															//UPDATE

        fin >> Z.Nsten;
        if (Z.Nsten > 0)
            Z.IDsten = new long [Z.Nsten];

        for ( long i = 0; i < Z.Nsten; i++ ){

                fin >> tmp_i;
                if (tmp_i == 1)
                    addBr = 0;																								//UPDATE																				
                else
                    addBr = Tr.NbrL + 0;																					//UPDATE

                if ((tmp_i != 1)&&(tmp_i != 2))
                    cout << "Wrong stenosis segment tree" << endl;

                fin >> tmp_i;
                Z.IDsten[i] = tmp_i + addBr;
                fin >> tmp_i;
                Tr.B[Z.IDsten[i] - 1].stenType = tmp_i; // stenType: 1 - stenosis, 2 - stent, 0 - usual vessel
        }



		//----- Making list of branches incoming and outgoing to every knot -----
		for ( long i = 0; i < Tr.Nkn; i++ ) {
			Tr.K [ i ].Nou = 0;
			Tr.K [ i ].Nin = 0;
		}
		long Iin, Iou;
		if ( ( Tr.ID == PULMVEN ) || ( Tr.ID == SYSVEN ) ) {
			//  Venous tree edges direction e1 <-- e2
			for ( long i = 0; i < Tr.Nbr; i++ ) {
				Tr.B [ i ].InvertPoints = 1;
				Iin = Tr.B [ i ].Kn1 -> ID;
				Iou = Tr.B [ i ].Kn2 -> ID;
				Tr.K [ Iin - 1 ].Nin += 1;
				Tr.K [ Iou - 1 ].Nou += 1;
				Tr.K [ Iin - 1 ].Bin.push_back( &Tr.B [ i ] );
				Tr.K [ Iou - 1 ].Bou.push_back( &Tr.B [ i ] );
			}
		} else {
			//  Arterial and Lymphatic trees edges direction e1 --> e2
			for ( long i = 0; i < Tr.Nbr; i++ ) {
				Tr.B [ i ].InvertPoints = 0;
				Iou = Tr.B [ i ].Kn1 -> ID;
				Iin = Tr.B [ i ].Kn2 -> ID;
				Tr.K [ Iin - 1 ].Nin += 1;
				Tr.K [ Iou - 1 ].Nou += 1;
				Tr.K [ Iin - 1 ].Bin.push_back( &Tr.B [ i ] );
				Tr.K [ Iou - 1 ].Bou.push_back( &Tr.B [ i ] );
			}
		}

		// Marking input and output as BIFURCATION
		for ( long i = 0; i < Tr.Nkn; i++ ) {
			if ( ( Tr.K [ i ].Nin + Tr.K [ i ].Nou ) == 1 ) {
				Tr.K [ i ].IG = FLOW;
			} else {
				Tr.K [ i ].IG = BIFURCATION;
			}
		}

        if (Z.Debug == 2)
        {
            cout << "Total knots loaded: " << Tr.Nkn << endl << "--------------------------------------" << endl;
            for ( long i = 0; i < Tr.Nkn; i++ ) printUzel ( Tr.K[i] );
            cout << "Total branches loaded: " << Tr.Nbr << endl << "--------------------------------------" << endl;
            for ( long i = 0; i < Tr.Nbr; i++ ) printVetv ( Tr.B[i] );
        }
        fin.close( );
	}

	inline void LoadMultiKnots( ) {  // not used in coronary
		long  KntID , szin , szou , LTrID;
		Derevo LTree;

		Globals::MKnots.multiKnotsfilename = trim( SharedDirectory ) + "multiknots.tre";
		ifstream fin( Globals::MKnots.multiKnotsfilename , ifstream::in );
		fin >> Globals::MKnots.sz;

		if ( Globals::MKnots.sz > 0 ) {
			Globals::MKnots.Lst.resize( Globals::MKnots.sz );
			for ( long i = 0; i < Globals::MKnots.sz; i++ ) {
				fin >> Globals::MKnots.Lst [ i ].sz;
				Globals::MKnots.Lst [ i ].GrpLst.resize( Globals::MKnots.Lst [ i ].sz );
				szin = 0;
				szou = 0;
				for ( long j = 0; j < Globals::MKnots.Lst [ i ].sz; j++ ) {
					fin >> MKnots.Lst [ i ].GrpLst [ j ].TrID;
					fin >> MKnots.Lst [ i ].GrpLst [ j ].sz;
					Globals::MKnots.Lst [ i ].GrpLst [ j ].KntLst.resize( Globals::MKnots.Lst [ i ].GrpLst [ j ].sz );
					for ( long k = 0; k < Globals::MKnots.Lst [ i ].GrpLst [ j ].sz; k++ ) {
						fin >> KntID;
						LTrID = MKnots.Lst [ i ].GrpLst [ j ].TrID;
						LTree = TreeLst [ LTrID ];
						Globals::MKnots.Lst [ i ].GrpLst [ j ].KntLst [ k ] = LTree.K [ KntID ];
						szin = szin + Globals::MKnots.Lst [ i ].GrpLst [ j ].KntLst [ k ].Nin;
						szou = szou + Globals::MKnots.Lst [ i ].GrpLst [ j ].KntLst [ k ].Nou;
					}
				}
				Globals::MKnots.Lst [ i ].szin = szin;
				Globals::MKnots.Lst [ i ].szou = szou;
			}
		}
		fin.close( );

		cout << "Total multiknots loaded: " << MKnots.sz << endl;
	}
};
