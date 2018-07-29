#define Rho1    2              //初期密度1
#define Vel1    0.5              //初期速度1
#define Press1  2.5              //初期圧力1
#define Rho2    1            //初期密度2
#define Vel2    -0.5            //初期速度2
#define Press2  2.5              //初期圧力2
#define EndTime 1.10           //終了時刻
#define XMax    1.0              //系の大きさ(xの最大値)
#define XMin    0.0             //系の大きさ(xの最小値)
#define YMax    1.022              //系の大きさ(yの最大値)
#define YMin    0.0             //系の大きさ(yの最小値)
#define DeltaX  (XMax - XMin) / 2  //系の大きさの半分
#define Xmid    (XMax + XMin) / 2   //系の中心座標
#define Alpha   1.0              //粘性項のα
#define Beta    2.0 * Alpha      //粘性項のβ
#define Epsilon 0.01             //粘性項のε
#define Neighbor_PTCL 4         //近傍粒子数
#define W0  0.025               //y方向速度の摂動
#define Lambda 1/2              //y方向速度の摂動
#define Sigma 0.05/sqrt(2)      //y方向速度の摂動
