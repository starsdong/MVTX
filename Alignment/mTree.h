typedef struct {
  Float_t ox;
  Float_t oy;
  Float_t oz;
  Float_t px;
  Float_t py;
  Float_t pz;
  Int_t   nh;
  Int_t   id[20]; // layer# [0-2] * 100 + sector#[0-12,16,20]
  Float_t xL_mc[20];
  Float_t yL_mc[20];
  Float_t zL_mc[20];
  Float_t xG_mc[20];
  Float_t yG_mc[20];
  Float_t zG_mc[20];
  Float_t xL_rc[20];
  Float_t yL_rc[20];
  Float_t zL_rc[20];
  Float_t xG_rc[20];
  Float_t yG_rc[20];
  Float_t zG_rc[20];  
} MVTXTREE;

