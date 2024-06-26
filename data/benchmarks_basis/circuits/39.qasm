OPENQASM 2.0;
include "qelib1.inc";
qreg q[10];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
cx q[1],q[4];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
cx q[4],q[1];
rz(-pi/4) q[1];
cx q[5],q[1];
rz(pi/4) q[1];
cx q[4],q[1];
rz(-pi/4) q[1];
rz(pi/4) q[4];
cx q[5],q[1];
rz(3*pi/4) q[1];
sx q[1];
rz(pi/2) q[1];
cx q[5],q[4];
rz(-pi/4) q[4];
rz(pi/4) q[5];
cx q[5],q[4];
cx q[1],q[4];
rz(-pi) q[1];
sx q[1];
rz(1.8551539769067755) q[1];
sx q[1];
rz(pi/2) q[6];
sx q[6];
rz(pi/2) q[6];
cx q[2],q[6];
rz(-pi/4) q[6];
cx q[2],q[6];
x q[2];
cx q[2],q[6];
rz(-pi/4) q[6];
cx q[2],q[6];
cx q[2],q[5];
cx q[5],q[2];
sx q[2];
rz(-pi/4) q[2];
rz(pi/2) q[5];
rz(pi/2) q[6];
sx q[6];
rz(pi/2) q[6];
rz(pi/4) q[7];
cx q[7],q[3];
rz(-pi/4) q[3];
cx q[7],q[3];
rz(3*pi/4) q[3];
sx q[3];
rz(pi/2) q[3];
cx q[8],q[0];
rz(1.99358308266723) q[0];
cx q[8],q[0];
rz(pi/2) q[8];
sx q[8];
rz(3*pi/4) q[8];
cx q[0],q[8];
rz(pi/4) q[8];
sx q[8];
rz(pi/2) q[8];
cx q[4],q[8];
rz(pi/4) q[8];
cx q[3],q[8];
rz(-pi/4) q[8];
cx q[4],q[8];
cx q[4],q[6];
rz(3.93663614920433) q[6];
cx q[4],q[6];
rz(pi/2) q[6];
sx q[6];
rz(pi/2) q[6];
rz(pi/4) q[8];
cx q[3],q[8];
cx q[3],q[1];
sx q[1];
rz(1.8551539769067746) q[1];
sx q[1];
rz(-pi) q[1];
cx q[3],q[1];
rz(3.1223865346471644) q[1];
sx q[1];
rz(-0.9305492109076567) q[1];
sx q[1];
rz(-0.5115483283636326) q[1];
rz(-pi/2) q[3];
cx q[4],q[3];
rz(pi/2) q[3];
cx q[2],q[3];
rz(pi/4) q[3];
cx q[2],q[3];
rz(-pi/4) q[3];
reset q[4];
rz(0.6406349125025548) q[4];
cx q[4],q[2];
rz(-0.6406349125025557) q[2];
cx q[4],q[2];
rz(0.6406349125025548) q[2];
rz(pi/2) q[4];
sx q[4];
rz(pi/2) q[4];
rz(pi/4) q[8];
sx q[8];
rz(3*pi/4) q[8];
cx q[0],q[8];
rz(1.31673914288684) q[0];
rz(-pi/4) q[8];
sx q[8];
rz(3*pi/4) q[8];
cx q[8],q[6];
rz(-pi/4) q[6];
cx q[8],q[6];
rz(pi/4) q[6];
sx q[6];
rz(pi/4) q[6];
rz(pi/2) q[9];
cx q[7],q[9];
cx q[9],q[7];
cx q[7],q[9];
rz(2.480658946512417) q[7];
sx q[7];
rz(-2.159782048177445) q[7];
cx q[7],q[5];
sx q[5];
rz(2.422866052018107) q[5];
sx q[5];
rz(-pi) q[5];
cx q[7],q[5];
sx q[5];
rz(-2.422866052018107) q[5];
sx q[5];
rz(-1.298169737933149) q[5];
rz(1.2798375987103) q[7];
cx q[8],q[7];
rz(-1.2798375987103) q[7];
cx q[8],q[7];
rz(-pi) q[7];
sx q[7];
rz(pi/2) q[7];
rz(pi/2) q[8];
sx q[8];
rz(pi/2) q[8];
cx q[3],q[8];
rz(-pi/4) q[8];
cx q[3],q[8];
x q[3];
cx q[3],q[8];
rz(-pi/4) q[8];
cx q[3],q[8];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
cx q[9],q[0];
rz(-1.31673914288684) q[0];
cx q[9],q[0];
rz(pi/2) q[9];
cx q[0],q[9];
sx q[9];
rz(2.367629145839784) q[9];
sx q[9];
rz(-pi) q[9];
cx q[0],q[9];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
cx q[0],q[5];
rz(-pi) q[0];
sx q[0];
rz(2.1590470371262453) q[0];
sx q[0];
sx q[5];
rz(2.1590470371262453) q[5];
sx q[5];
rz(-pi) q[5];
cx q[0],q[5];
rz(pi/2) q[0];
sx q[0];
rz(2.9933825361614588) q[0];
cx q[1],q[0];
rz(0.49234220942100615) q[0];
sx q[0];
rz(-0.20018248489235546) q[0];
sx q[0];
cx q[1],q[0];
sx q[0];
rz(-0.20018248489235546) q[0];
sx q[0];
rz(2.0120623981996726) q[0];
cx q[0],q[4];
cx q[1],q[8];
rz(-pi/4) q[4];
cx q[0],q[4];
rz(pi/2) q[0];
cx q[1],q[0];
sx q[0];
rz(2.800795836167289) q[0];
sx q[0];
rz(-pi) q[0];
cx q[1],q[0];
rz(-3.078822409764104) q[0];
sx q[0];
rz(-1.9109726142893377) q[0];
sx q[0];
rz(-2.5171073092499654) q[0];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3*pi/4) q[4];
sx q[4];
rz(pi/2) q[4];
rz(-3.1351194857452915) q[5];
rz(-pi/2) q[8];
sx q[8];
rz(pi/2) q[8];
sx q[9];
rz(-2.367629145839784) q[9];
sx q[9];
rz(pi/2) q[9];
cx q[9],q[6];
rz(pi/4) q[6];
sx q[6];
rz(-0.7805873884246095) q[6];
cx q[6],q[5];
rz(-2.3610052651651854) q[5];
sx q[5];
rz(-3.010700358035768) q[5];
sx q[5];
cx q[6],q[5];
sx q[5];
rz(-3.010700358035768) q[5];
sx q[5];
rz(2.0819055084589344) q[5];
cx q[5],q[3];
rz(-pi/4) q[3];
cx q[5],q[3];
x q[5];
cx q[5],q[3];
rz(-pi/4) q[3];
cx q[5],q[3];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
x q[5];
rz(-0.2743886039980423) q[6];
sx q[6];
rz(-1.2267720821494166) q[6];
sx q[6];
rz(2.9247877272699645) q[6];
cx q[0],q[6];
rz(-pi/4) q[6];
cx q[0],q[6];
x q[0];
cx q[0],q[6];
rz(-pi/4) q[6];
cx q[0],q[6];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[6];
sx q[6];
rz(pi/2) q[6];
rz(pi/2) q[9];
cx q[7],q[9];
cx q[9],q[7];
rz(2.93464255530042) q[7];
cx q[2],q[7];
rz(-2.93464255530042) q[7];
cx q[2],q[7];
rz(pi/2) q[7];
sx q[7];
rz(pi/2) q[7];
cx q[3],q[7];
rz(-pi/4) q[7];
cx q[4],q[7];
rz(pi/4) q[7];
cx q[3],q[7];
rz(pi/4) q[3];
rz(-pi/4) q[7];
cx q[4],q[7];
cx q[4],q[3];
rz(-pi/4) q[3];
rz(pi/4) q[4];
cx q[4],q[3];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
cx q[4],q[5];
cx q[5],q[4];
rz(pi/2) q[4];
sx q[4];
rz(-pi/2) q[4];
rz(2.55525691328393) q[5];
sx q[5];
rz(7.53193533557614) q[5];
sx q[5];
rz(13.215299693318) q[5];
cx q[5],q[4];
rz(-pi) q[4];
sx q[4];
rz(pi/2) q[4];
rz(3.020897015567341) q[7];
sx q[7];
rz(-1.1984429840240196) q[7];
sx q[7];
rz(-1.014280253205305) q[7];
cx q[8],q[3];
rz(-pi/4) q[3];
cx q[8],q[3];
x q[8];
cx q[8],q[3];
rz(-pi/4) q[3];
cx q[8],q[3];
rz(pi/4) q[8];
cx q[8],q[3];
rz(-pi/4) q[3];
cx q[8],q[3];
rz(pi/4) q[3];
sx q[3];
rz(pi/4) q[3];
cx q[6],q[3];
rz(pi/4) q[3];
sx q[3];
rz(1.3935958638244097) q[6];
sx q[6];
rz(pi/2) q[6];
cx q[7],q[6];
rz(pi/2) q[6];
sx q[6];
rz(pi/2) q[6];
rz(2.401075441541785) q[7];
rz(pi/2) q[8];
sx q[8];
rz(pi/2) q[8];
cx q[0],q[8];
rz(0.646358515983621) q[8];
cx q[0],q[8];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/4) q[8];
cx q[5],q[8];
rz(pi/4) q[8];
sx q[8];
rz(pi/2) q[8];
cx q[3],q[8];
rz(pi/4) q[8];
cx q[0],q[8];
rz(-pi/4) q[8];
cx q[3],q[8];
rz(pi/4) q[8];
cx q[0],q[8];
rz(pi/4) q[8];
sx q[8];
rz(3*pi/4) q[8];
cx q[5],q[8];
cx q[5],q[3];
cx q[3],q[5];
rz(-3*pi/4) q[5];
sx q[5];
rz(0.4829104462353926) q[5];
sx q[5];
cx q[7],q[5];
sx q[5];
rz(0.4829104462353917) q[5];
sx q[5];
rz(-pi) q[5];
cx q[7],q[5];
rz(-2.4646132807455334) q[5];
rz(pi/4) q[8];
sx q[8];
rz(pi/2) q[8];
rz(pi/2) q[9];
sx q[9];
rz(3*pi/4) q[9];
cx q[2],q[9];
cx q[9],q[2];
rz(-1.8914993994962366) q[2];
sx q[2];
rz(1.1615217699895606) q[2];
sx q[2];
cx q[4],q[2];
rz(1.9428253816608) q[2];
cx q[4],q[2];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
cx q[2],q[8];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
rz(-pi) q[4];
sx q[4];
rz(-pi) q[4];
cx q[0],q[4];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
cx q[4],q[0];
rz(-pi/4) q[0];
cx q[8],q[2];
rz(-pi/4) q[2];
cx q[9],q[1];
rz(1.17594775294878) q[1];
cx q[9],q[1];
rz(pi/2) q[1];
sx q[1];
rz(-2.945088085517929) q[1];
rz(-pi/4) q[9];
sx q[9];
rz(pi/2) q[9];
cx q[9],q[1];
sx q[1];
rz(1.7287298519797103) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[9];
sx q[9];
rz(1.7287298519797103) q[9];
sx q[9];
cx q[9],q[1];
rz(-0.9819027314693152) q[1];
cx q[1],q[0];
rz(pi/4) q[0];
cx q[4],q[0];
rz(-pi/4) q[0];
cx q[1],q[0];
rz(3*pi/4) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/4) q[4];
cx q[1],q[4];
rz(pi/4) q[1];
rz(-pi/4) q[4];
cx q[1],q[4];
cx q[0],q[4];
rz(-1.415054261973678) q[0];
sx q[0];
rz(-1.3154962582966778) q[0];
sx q[0];
rz(3*pi/4) q[0];
sx q[1];
rz(pi/2) q[4];
sx q[4];
rz(pi/2) q[4];
cx q[6],q[4];
rz(-pi/4) q[4];
cx q[3],q[4];
rz(pi/4) q[4];
cx q[6],q[4];
rz(-pi/4) q[4];
cx q[3],q[4];
rz(3*pi/4) q[4];
sx q[4];
rz(pi/2) q[4];
rz(pi/4) q[6];
cx q[3],q[6];
rz(pi/4) q[3];
rz(-pi/4) q[6];
cx q[3],q[6];
rz(-2.650448616720265) q[6];
sx q[6];
rz(-2.2259808941140538) q[6];
sx q[6];
rz(1.655046421146265) q[6];
rz(pi/2) q[9];
sx q[9];
rz(-pi/2) q[9];
cx q[9],q[2];
rz(pi/4) q[2];
cx q[8],q[2];
rz(-pi/4) q[2];
rz(pi/4) q[8];
cx q[9],q[2];
rz(3*pi/4) q[2];
sx q[2];
rz(pi/2) q[2];
cx q[9],q[8];
rz(-pi/4) q[8];
rz(pi/4) q[9];
cx q[9],q[8];
cx q[2],q[8];
rz(2.34449073805889) q[2];
sx q[2];
rz(3.98563385895264) q[2];
sx q[2];
rz(12.3446157203284) q[2];
cx q[2],q[3];
rz(1.31325172667429) q[3];
cx q[2],q[3];
rz(-2.748602092011037) q[3];
cx q[3],q[7];
rz(-pi/4) q[7];
cx q[3],q[7];
rz(2.1788945584915522) q[3];
rz(-0.11977913445605726) q[7];
sx q[8];
cx q[8],q[1];
rz(2.94956184201571) q[1];
cx q[8],q[1];
sx q[1];
rz(0.922085367209923) q[1];
cx q[5],q[1];
rz(-0.6769793728442606) q[1];
sx q[1];
rz(-2.5031315634186457) q[1];
sx q[1];
cx q[5],q[1];
sx q[1];
rz(-2.5031315634186457) q[1];
sx q[1];
rz(0.04435608584627815) q[1];
cx q[3],q[1];
rz(1.7709901698607649) q[1];
sx q[1];
rz(-0.9900067693953396) q[1];
sx q[1];
cx q[3],q[1];
sx q[1];
rz(-0.9900067693953396) q[1];
sx q[1];
rz(-0.48965592327781193) q[1];
rz(-pi) q[3];
sx q[3];
rz(3*pi/4) q[3];
cx q[5],q[0];
rz(pi/4) q[0];
sx q[0];
rz(-0.01928216477745126) q[0];
rz(pi/2) q[5];
sx q[5];
rz(pi/2) q[5];
rz(-pi) q[8];
sx q[8];
rz(-pi) q[8];
x q[9];
cx q[4],q[9];
rz(-pi/4) q[9];
cx q[8],q[9];
rz(pi/4) q[9];
cx q[4],q[9];
rz(pi/4) q[4];
rz(-pi/4) q[9];
cx q[8],q[9];
cx q[8],q[4];
rz(-pi/4) q[4];
rz(pi/4) q[8];
cx q[8],q[4];
cx q[4],q[8];
rz(pi/2) q[4];
sx q[4];
rz(pi/2) q[4];
cx q[8],q[4];
rz(-pi/4) q[4];
cx q[2],q[4];
rz(pi/4) q[4];
cx q[8],q[4];
rz(-pi/4) q[4];
cx q[2],q[4];
rz(3*pi/4) q[4];
sx q[4];
rz(pi/2) q[4];
rz(pi/4) q[8];
cx q[2],q[8];
rz(pi/4) q[2];
rz(-pi/4) q[8];
cx q[2],q[8];
rz(pi/2) q[2];
cx q[4],q[8];
cx q[4],q[2];
sx q[2];
rz(0.9985401411523913) q[2];
sx q[2];
rz(-pi) q[2];
cx q[4],q[2];
rz(2.375967625155635) q[2];
sx q[2];
rz(-0.7292030811809997) q[2];
sx q[2];
rz(-3.0476395150483215) q[2];
cx q[2],q[3];
rz(3*pi/4) q[3];
sx q[3];
rz(-pi) q[3];
rz(0.159443700725209) q[4];
sx q[4];
rz(6.4276477978535) q[4];
sx q[4];
rz(9.26533426004417) q[4];
cx q[4],q[1];
sx q[1];
rz(2.4132991927479033) q[1];
sx q[1];
rz(-pi) q[1];
cx q[4],q[1];
sx q[1];
rz(-2.4132991927479033) q[1];
sx q[1];
rz(pi/2) q[1];
cx q[1],q[3];
rz(0.76071178026888) q[3];
cx q[1],q[3];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
x q[4];
rz(pi/2) q[8];
sx q[8];
rz(pi/2) q[8];
rz(-pi/4) q[9];
sx q[9];
rz(pi/2) q[9];
cx q[9],q[8];
rz(0.213569208374314) q[8];
cx q[9],q[8];
rz(pi/2) q[8];
sx q[8];
rz(1.0005558335079074) q[8];
cx q[7],q[8];
rz(0.9807206108432638) q[8];
sx q[8];
rz(-2.4280059538784213) q[8];
sx q[8];
cx q[7],q[8];
cx q[7],q[5];
rz(pi/2) q[5];
sx q[5];
rz(pi/2) q[5];
rz(-pi) q[7];
sx q[7];
rz(1.6358562705791737) q[7];
sx q[7];
cx q[2],q[7];
sx q[7];
rz(1.6358562705791737) q[7];
sx q[7];
rz(-pi) q[7];
cx q[2],q[7];
rz(2.43550676524948) q[2];
sx q[2];
rz(pi/2) q[2];
rz(0.0005916203305114109) q[7];
sx q[8];
rz(-2.4280059538784213) q[8];
sx q[8];
rz(-0.3053422734019473) q[8];
cx q[8],q[0];
rz(3.0364548094354635) q[0];
sx q[0];
rz(-1.0821796754035056) q[0];
sx q[0];
cx q[8],q[0];
sx q[0];
rz(-1.0821796754035056) q[0];
sx q[0];
rz(-3.0171726446580127) q[0];
cx q[0],q[6];
sx q[6];
rz(2.9740853044030944) q[6];
sx q[6];
rz(-pi) q[6];
cx q[0],q[6];
rz(-0.8213112892589471) q[0];
cx q[4],q[0];
sx q[0];
rz(1.6189117613822939) q[0];
sx q[0];
rz(-pi) q[0];
sx q[4];
rz(1.6189117613822939) q[4];
sx q[4];
rz(-pi) q[4];
cx q[4],q[0];
rz(-0.602285446686444) q[0];
cx q[2],q[0];
sx q[0];
rz(1.0014269939824825) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(1.0014269939824825) q[2];
sx q[2];
rz(-pi) q[2];
cx q[2],q[0];
rz(2.002576484500681) q[0];
sx q[0];
rz(-1.1913775418933916) q[0];
sx q[0];
rz(-2.1375696628815337) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
rz(pi/2) q[4];
sx q[4];
rz(-pi/2) q[4];
rz(1.07904391606281) q[8];
sx q[8];
rz(8.86952086795193) q[8];
sx q[8];
rz(13.0223452177217) q[8];
reset q[8];
rz(5.07564626791568) q[8];
sx q[8];
rz(7.17354284334409) q[8];
sx q[8];
rz(15.6054223686819) q[8];
rz(-0.34489490022391234) q[9];
sx q[9];
rz(3*pi/4) q[9];
cx q[9],q[5];
rz(-pi/4) q[5];
cx q[9],q[5];
rz(pi/4) q[5];
cx q[5],q[6];
rz(pi/4) q[5];
cx q[5],q[7];
x q[6];
rz(-pi/4) q[7];
cx q[5],q[7];
rz(pi/2) q[5];
rz(3*pi/4) q[7];
sx q[7];
rz(3*pi/4) q[7];
cx q[8],q[6];
cx q[6],q[8];
rz(1.6686828984697737) q[6];
sx q[6];
rz(-0.12091316100130678) q[6];
rz(pi/2) q[8];
sx q[8];
rz(pi/2) q[8];
sx q[9];
cx q[9],q[1];
rz(3.95713558116232) q[1];
cx q[9],q[1];
rz(-pi) q[1];
sx q[1];
rz(-1.518346628445467) q[1];
cx q[1],q[4];
rz(-1.6232460251443275) q[4];
cx q[1],q[4];
cx q[1],q[7];
rz(-1.518346628445467) q[4];
sx q[4];
rz(pi/2) q[4];
cx q[4],q[5];
cx q[5],q[4];
rz(pi/2) q[4];
sx q[4];
rz(pi/2) q[4];
cx q[0],q[4];
rz(-pi/4) q[4];
cx q[0],q[4];
rz(pi/2) q[0];
sx q[0];
rz(4.79453387945117) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(3*pi/4) q[4];
sx q[4];
rz(pi/2) q[4];
rz(-pi/4) q[7];
cx q[2],q[7];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
cx q[5],q[2];
rz(3.55553068457032) q[2];
cx q[5],q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/4) q[2];
sx q[5];
rz(2.2645711871646537) q[5];
sx q[5];
rz(pi/2) q[5];
rz(pi/4) q[7];
cx q[1],q[7];
cx q[1],q[8];
rz(pi/4) q[7];
sx q[7];
rz(pi/2) q[7];
rz(-pi/4) q[8];
cx q[1],q[8];
x q[1];
cx q[1],q[8];
rz(-pi/4) q[8];
cx q[1],q[8];
x q[1];
rz(-pi/2) q[8];
sx q[8];
rz(-pi/2) q[8];
rz(-2.970534243820019) q[9];
sx q[9];
rz(-0.6165372548940695) q[9];
sx q[9];
rz(-1.2420911562626014) q[9];
cx q[3],q[9];
rz(-pi/4) q[9];
cx q[3],q[9];
x q[3];
cx q[3],q[9];
rz(-pi/4) q[9];
cx q[3],q[9];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
cx q[7],q[3];
rz(-pi/4) q[3];
cx q[7],q[3];
x q[7];
cx q[7],q[3];
rz(-pi/4) q[3];
cx q[7],q[3];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
cx q[3],q[0];
rz(3.81715645335531) q[0];
cx q[3],q[0];
rz(pi/2) q[7];
sx q[7];
rz(pi/2) q[7];
cx q[4],q[7];
reset q[4];
rz(pi/4) q[4];
cx q[8],q[7];
rz(-pi/4) q[7];
cx q[6],q[7];
rz(pi/4) q[7];
cx q[8],q[7];
rz(-pi/4) q[7];
cx q[6],q[7];
rz(pi/4) q[7];
cx q[3],q[7];
rz(-pi/4) q[7];
cx q[0],q[7];
rz(pi/4) q[7];
cx q[3],q[7];
rz(pi/4) q[3];
rz(-pi/4) q[7];
cx q[0],q[7];
cx q[0],q[3];
rz(pi/4) q[0];
rz(-pi/4) q[3];
cx q[0],q[3];
rz(pi/2) q[3];
sx q[3];
rz(3*pi/4) q[3];
rz(1.9012125203978183) q[7];
sx q[7];
rz(-0.5745408728373533) q[7];
sx q[7];
rz(1.3181573960589033) q[7];
rz(pi/4) q[8];
cx q[6],q[8];
rz(pi/4) q[6];
rz(-pi/4) q[8];
cx q[6],q[8];
sx q[6];
cx q[1],q[6];
rz(0.716738123016963) q[6];
cx q[1],q[6];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[6];
rz(-pi) q[6];
cx q[6],q[7];
rz(2.93764621600944) q[7];
cx q[6],q[7];
rz(-0.4687902723741755) q[7];
rz(pi/2) q[8];
rz(-pi/2) q[9];
sx q[9];
rz(-0.8030818493604173) q[9];
sx q[9];
rz(-pi/2) q[9];
cx q[2],q[9];
rz(-pi/4) q[9];
cx q[2],q[9];
rz(3.1015495754427356) q[2];
cx q[2],q[5];
rz(-3.101549575442736) q[5];
cx q[2],q[5];
cx q[4],q[2];
rz(-pi/4) q[2];
cx q[4],q[2];
rz(pi/4) q[2];
rz(3.1015495754427356) q[5];
cx q[5],q[3];
rz(-pi/4) q[3];
cx q[0],q[3];
rz(2.848889717287685) q[0];
cx q[0],q[2];
rz(-2.848889717287686) q[2];
cx q[0],q[2];
rz(1.6195240066697583) q[0];
rz(-1.8634992630970046) q[2];
sx q[2];
rz(pi/2) q[2];
rz(pi/4) q[3];
cx q[5],q[3];
rz(-1.4001040279463286) q[3];
sx q[3];
rz(-1.7467553335000208) q[3];
sx q[3];
rz(0.5268467235168366) q[3];
rz(pi/2) q[5];
sx q[5];
rz(3*pi/4) q[5];
cx q[7],q[0];
rz(0.46311177419375893) q[0];
sx q[0];
rz(-0.1521953102605842) q[0];
sx q[0];
cx q[7],q[0];
sx q[0];
rz(-0.1521953102605842) q[0];
sx q[0];
rz(-2.0826357808635176) q[0];
rz(-pi/4) q[9];
sx q[9];
rz(-3*pi/4) q[9];
sx q[9];
rz(-pi/4) q[9];
cx q[8],q[9];
rz(-pi/4) q[9];
cx q[1],q[9];
cx q[1],q[5];
rz(-pi/4) q[5];
rz(pi/4) q[9];
cx q[8],q[9];
rz(2.5822280228300016) q[8];
cx q[8],q[4];
rz(-2.582228022830003) q[4];
cx q[8],q[4];
rz(2.5822280228300016) q[4];
cx q[6],q[4];
rz(pi/2) q[6];
sx q[6];
rz(pi/2) q[6];
cx q[4],q[6];
rz(-pi/4) q[6];
cx q[8],q[6];
rz(pi/4) q[6];
cx q[4],q[6];
rz(pi/4) q[4];
rz(-pi/4) q[6];
cx q[8],q[6];
rz(3*pi/4) q[6];
sx q[6];
rz(pi/2) q[6];
cx q[8],q[4];
rz(-pi/4) q[4];
rz(pi/4) q[8];
cx q[8],q[4];
cx q[6],q[4];
rz(pi/4) q[9];
sx q[9];
rz(pi/2) q[9];
cx q[9],q[5];
rz(pi/4) q[5];
cx q[1],q[5];
sx q[1];
rz(pi/4) q[5];
sx q[5];
rz(-pi) q[5];
cx q[9],q[2];
rz(5.75450394975605) q[2];
cx q[9],q[2];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
