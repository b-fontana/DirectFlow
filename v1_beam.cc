

#include "/Users/marialex/STAR/Utils/functions.h"
#include "./functions_proton_spec.h"



void v1_beam()
{
    cout << "v1_beam macro started" << endl;
    gRandom->SetSeed(0);
    gStyle->SetPalette(56); // 53 = black body radiation, 56 = inverted black body radiator, 103 = sunset, 87 == light temperature


    TVector3 vec_p;
    vec_p.SetXYZ(0.0,0.0,5000.0); // 5000
    Double_t pos_xyz[3] = {0.0,-0.0,-6990.0};

    TLorentzVector* TLV_particle;
    TLV_particle = new TLorentzVector();

    Double_t particle_energy = TMath::Sqrt(vec_p.Mag()*vec_p.Mag() + 0.938*0.938);
    TLV_particle       ->SetPxPyPzE(vec_p.x(),vec_p.y(),vec_p.z(),particle_energy); // nucleon momentum in cm frame, Fermi momentum + v1 momentum

    TEveManager::Create();
    TEveLine* TEveLine_beam_axis = NULL;
    TEveLine_beam_axis = new TEveLine();
    TEveLine_beam_axis ->SetNextPoint(0.0,0.0,-7650.0);
    TEveLine_beam_axis ->SetNextPoint(0.0,0.0,7650.0);
    TEveLine_beam_axis ->SetName("beam axis");
    TEveLine_beam_axis ->SetLineStyle(1);
    TEveLine_beam_axis ->SetLineWidth(4);
    TEveLine_beam_axis ->SetMainAlpha(0.7);
    TEveLine_beam_axis ->SetMainColor(kBlue);
    gEve->AddElement(TEveLine_beam_axis);

    vector<TEveBox*> vec_eve_magnets;
    vec_eve_magnets.resize(12);
    TString label_magnet_names[12] = {"D1_neg","Q4_neg","Q3_neg","Q2_neg","Q1_neg","D_corr","Muon","Q1_pos","Q2_pos","Q3_pos","Q4_pos","D1_pos"};
    Double_t z_start_magnets[12] = {-5840.0 - 945.0,-4730.0 - 630.0,-3830.0 - 550.0,-3180.0 - 550.0,-2300.0 - 630.0,-1920.0 - 190.0,-750.0 - 430.0,
    2300.0,3180.0,3830.0,4730.0,5840.0};
    Double_t z_stop_magnets[12] = {-5840.0,-4730.0,-3830.0,-3180.0,-2300.0,-1920.0,-750.0,
    2300.0 + 630.0,3180.0 + 550.0,3830.0 + 550.0,4730.0 + 630.0,5840.0 + 945.0};
    Double_t x_size_magnets = 100.0;
    Double_t y_size_magnets = 100.0;
    Int_t color_magnets[12] = {kBlue,kYellow,kYellow,kYellow,kYellow,kBlue+1,kMagenta,kYellow,kYellow,kYellow,kYellow,kBlue};
    for(Int_t i_magnet = 0; i_magnet < 12; i_magnet++)
    {
        vec_eve_magnets[i_magnet] = new TEveBox;
        vec_eve_magnets[i_magnet] ->SetName(label_magnet_names[i_magnet].Data());

        vec_eve_magnets[i_magnet] ->SetVertex(0,-x_size_magnets,-y_size_magnets,z_start_magnets[i_magnet]);
        vec_eve_magnets[i_magnet] ->SetVertex(1,x_size_magnets,-y_size_magnets,z_start_magnets[i_magnet]);
        vec_eve_magnets[i_magnet] ->SetVertex(2,x_size_magnets,y_size_magnets,z_start_magnets[i_magnet]);
        vec_eve_magnets[i_magnet] ->SetVertex(3,-x_size_magnets,y_size_magnets,z_start_magnets[i_magnet]);
        vec_eve_magnets[i_magnet] ->SetVertex(4,-x_size_magnets,-y_size_magnets,z_stop_magnets[i_magnet]);
        vec_eve_magnets[i_magnet] ->SetVertex(5,x_size_magnets,-y_size_magnets,z_stop_magnets[i_magnet]);
        vec_eve_magnets[i_magnet] ->SetVertex(6,x_size_magnets,y_size_magnets,z_stop_magnets[i_magnet]);
        vec_eve_magnets[i_magnet] ->SetVertex(7,-x_size_magnets,y_size_magnets,z_stop_magnets[i_magnet]);

        vec_eve_magnets[i_magnet]->SetMainColor(color_magnets[i_magnet]);
        vec_eve_magnets[i_magnet]->SetMainTransparency(75); // the higher the value the more transparent
        gEve->AddElement(vec_eve_magnets[i_magnet]);
    }

    vector<TEveLine*> vec_particle_track;
    vec_particle_track.resize(1);
    vec_particle_track[0] = new TEveLine();



    //-------------------------------------------------------------------------------------------
    // Scanned sigma x (m) as a function of z (m) of the LHC beam  (John Jowett)
    Double_t z_LHC[37]       = {-67.643,-58.738,-54.063,-50.500,-47.161,-44.267,-40.259,-38.701,-32.244,-30.241,-26.679,-23.784,-20.222,-17.105,-13.098,-9.3135,-4.6382,-0.4081,3.37662,8.05195,13.6178,17.6252,20.2968,22.7458,26.0853,28.9796,31.8738,34.7681,38.1076,40.7792,44.1187,46.3451,49.2393,53.0241,58.3673,63.4879,69.0538};
    Double_t sigma_x_LHC[37] = {0.0013124,0.0014062,0.0014562,0.0013687,0.0012687,0.0010499,0.0008562,0.0007999,0.0007562,0.0007812,0.0007624,0.0007437,0.0006249,0.0005249,0.0004062,0.0002812,0.0001437,4.37485e-05,0.0001312,0.0002687,0.0004374,0.0005624,0.0006562,0.0007187,0.0009124,0.0010812,0.0013124,0.0014312,0.0015624,0.0014874,0.0014249,0.0012374,0.0011374,0.0010187,0.0010062,0.0009812,0.0009687};
    //-------------------------------------------------------------------------------------------



    Int_t N_steps_used = 0;
    Track_eta_pT_phi_q_xyz_m(TLV_particle->Px(),TLV_particle->Py(),TLV_particle->Pz(),0,pos_xyz[0],pos_xyz[1],pos_xyz[2],0.938,
                             N_total_steps,step_size,N_steps_used); // fills track_pos_dir

    for(Int_t i_step = 0; i_step < N_steps_used; i_step++)
    {
        vec_particle_track[0]        ->SetNextPoint(track_pos_dir[i_step][0],track_pos_dir[i_step][1],track_pos_dir[i_step][2]);
        //printf("i_step: %d, pos: {%4.3f, %4.3f, %4.3f} \n",i_step,track_pos_dir[i_step][0],track_pos_dir[i_step][1],track_pos_dir[i_step][2]);

    }

    HistName = "track ";
    HistName += 0;
    vec_particle_track[0] ->SetName(HistName.Data());
    vec_particle_track[0] ->SetLineStyle(1);
    vec_particle_track[0] ->SetLineWidth(5);
    vec_particle_track[0] ->SetMainColor(kRed);
    vec_particle_track[0] ->SetMainAlpha(0.7);
    gEve->AddElement(vec_particle_track[0]);




    gEve->Redraw3D(kTRUE);
}

