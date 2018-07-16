#define Rho1    1.0              //初期密度1
#define Vel1    0              //初期速度1
#define Press1  1.0              //初期圧力1
#define Rho2    0.125            //初期密度2
#define Vel2    0             //初期速度2
#define Press2  0.1              //初期圧力2
#define EndTime 0.25            //終了時刻
#define XMax    1.0              //系の大きさ(xの最大値)
#define XMin    0             //系の大きさ(xの最小値)
#define DeltaX  (XMax - XMin) / 2  //系の大きさの半分
#define Xmid    (XMax + XMin) / 2   //系の中心座標
#define Alpha   1.0              //粘性項のα
#define Beta    2.0 * Alpha      //粘性項のβ
#define Epsilon 0.01             //粘性項のε
#define Neighbor_PTCL 1          //近傍粒子数
