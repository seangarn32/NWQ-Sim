OPENQASM 2.0;
include "qelib1.inc";
qreg q37[24];
creg c2[24];
x q37[0];
x q37[2];
x q37[4];
x q37[6];
x q37[8];
x q37[10];
barrier q37[0],q37[1],q37[2],q37[3],q37[4],q37[5],q37[6],q37[7],q37[8],q37[9],q37[10],q37[11],q37[12],q37[13],q37[14],q37[15],q37[16],q37[17],q37[18],q37[19],q37[20],q37[21],q37[22],q37[23];
sx q37[0];
u2(-pi,-pi/2) q37[1];
cx q37[1],q37[0];
rz(-0.9853981633974493) q37[0];
u2(0.9853981633974493,pi/2) q37[1];
cx q37[1],q37[0];
rz(0.9853981633974496) q37[0];
sx q37[1];
cx q37[1],q37[0];
u2(-pi/2,-pi/2) q37[0];
u2(pi/2,0) q37[1];
u2(-pi,-pi/2) q37[10];
sx q37[11];
u2(-pi,-pi/2) q37[14];
u2(-pi,-pi/2) q37[15];
u2(-pi,-pi/2) q37[18];
u2(-pi,-pi/2) q37[19];
u2(-pi,-pi/2) q37[2];
cx q37[2],q37[1];
rz(-0.8000000000000006) q37[1];
u2(0.8000000000000007,pi/2) q37[2];
cx q37[2],q37[1];
rz(0.8000000000000008) q37[1];
sx q37[2];
cx q37[2],q37[1];
u2(0,pi/2) q37[1];
cx q37[1],q37[0];
rz(0.5853981633974487) q37[0];
u2(0.5853981633974485,pi/2) q37[1];
cx q37[1],q37[0];
rz(0.5853981633974479) q37[0];
sx q37[1];
cx q37[1],q37[0];
rz(-pi) q37[0];
u3(pi,-3*pi/4,3*pi/4) q37[1];
u3(pi,-3*pi/4,3*pi/4) q37[2];
u2(-pi,-pi/2) q37[22];
cx q37[22],q37[11];
rz(-0.9853981633974493) q37[11];
u2(0.9853981633974493,pi/2) q37[22];
cx q37[22],q37[11];
rz(0.9853981633974496) q37[11];
sx q37[22];
cx q37[22],q37[11];
u2(-pi/2,-pi/2) q37[11];
u2(pi/2,0) q37[22];
u2(-pi,-pi/2) q37[23];
cx q37[23],q37[22];
rz(-0.8000000000000006) q37[22];
u2(0.8000000000000007,pi/2) q37[23];
cx q37[23],q37[22];
rz(0.8000000000000008) q37[22];
sx q37[23];
cx q37[23],q37[22];
u2(0,pi/2) q37[22];
cx q37[22],q37[11];
rz(0.5853981633974487) q37[11];
u2(0.5853981633974485,pi/2) q37[22];
cx q37[22],q37[11];
rz(0.5853981633974479) q37[11];
sx q37[22];
cx q37[22],q37[11];
rz(-pi) q37[11];
u3(pi,-3*pi/4,3*pi/4) q37[22];
u3(pi,-3*pi/4,3*pi/4) q37[23];
sx q37[3];
cx q37[14],q37[3];
u2(0.9853981633974493,pi/2) q37[14];
rz(-0.9853981633974493) q37[3];
cx q37[14],q37[3];
sx q37[14];
rz(0.9853981633974496) q37[3];
cx q37[14],q37[3];
u2(pi/2,0) q37[14];
cx q37[15],q37[14];
rz(-0.8000000000000006) q37[14];
u2(0.8000000000000007,pi/2) q37[15];
cx q37[15],q37[14];
rz(0.8000000000000008) q37[14];
sx q37[15];
cx q37[15],q37[14];
u2(0,pi/2) q37[14];
u3(pi,-3*pi/4,3*pi/4) q37[15];
u2(-pi/2,-pi/2) q37[3];
cx q37[14],q37[3];
u2(0.5853981633974485,pi/2) q37[14];
rz(0.5853981633974487) q37[3];
cx q37[14],q37[3];
sx q37[14];
rz(0.5853981633974479) q37[3];
cx q37[14],q37[3];
u3(pi,-3*pi/4,3*pi/4) q37[14];
rz(-pi) q37[3];
sx q37[4];
u2(-pi,-pi/2) q37[5];
cx q37[5],q37[4];
rz(-0.9853981633974493) q37[4];
u2(0.9853981633974493,pi/2) q37[5];
cx q37[5],q37[4];
rz(0.9853981633974496) q37[4];
sx q37[5];
cx q37[5],q37[4];
u2(-pi/2,-pi/2) q37[4];
u2(pi/2,0) q37[5];
u2(-pi,-pi/2) q37[6];
cx q37[6],q37[5];
rz(-0.8000000000000006) q37[5];
u2(0.8000000000000007,pi/2) q37[6];
cx q37[6],q37[5];
rz(0.8000000000000008) q37[5];
sx q37[6];
cx q37[6],q37[5];
u2(0,pi/2) q37[5];
cx q37[5],q37[4];
rz(0.5853981633974487) q37[4];
u2(0.5853981633974485,pi/2) q37[5];
cx q37[5],q37[4];
rz(0.5853981633974479) q37[4];
sx q37[5];
cx q37[5],q37[4];
rz(-pi) q37[4];
u3(pi,-3*pi/4,3*pi/4) q37[5];
u3(pi,-3*pi/4,3*pi/4) q37[6];
sx q37[7];
cx q37[18],q37[7];
u2(0.9853981633974493,pi/2) q37[18];
rz(-0.9853981633974493) q37[7];
cx q37[18],q37[7];
sx q37[18];
rz(0.9853981633974496) q37[7];
cx q37[18],q37[7];
u2(pi/2,0) q37[18];
cx q37[19],q37[18];
rz(-0.8000000000000006) q37[18];
u2(0.8000000000000007,pi/2) q37[19];
cx q37[19],q37[18];
rz(0.8000000000000008) q37[18];
sx q37[19];
cx q37[19],q37[18];
u2(0,pi/2) q37[18];
u3(pi,-3*pi/4,3*pi/4) q37[19];
u2(-pi/2,-pi/2) q37[7];
cx q37[18],q37[7];
u2(0.5853981633974485,pi/2) q37[18];
rz(0.5853981633974487) q37[7];
cx q37[18],q37[7];
sx q37[18];
rz(0.5853981633974479) q37[7];
cx q37[18],q37[7];
u3(pi,-3*pi/4,3*pi/4) q37[18];
rz(-pi) q37[7];
sx q37[8];
u2(-pi,-pi/2) q37[9];
cx q37[9],q37[8];
rz(-0.9853981633974493) q37[8];
u2(0.9853981633974493,pi/2) q37[9];
cx q37[9],q37[8];
rz(0.9853981633974496) q37[8];
sx q37[9];
cx q37[9],q37[8];
u2(-pi/2,-pi/2) q37[8];
u2(pi/2,0) q37[9];
cx q37[10],q37[9];
u2(0.8000000000000007,pi/2) q37[10];
rz(-0.8000000000000006) q37[9];
cx q37[10],q37[9];
sx q37[10];
rz(0.8000000000000008) q37[9];
cx q37[10],q37[9];
u3(pi,-3*pi/4,3*pi/4) q37[10];
u2(0,pi/2) q37[9];
cx q37[9],q37[8];
rz(0.5853981633974487) q37[8];
u2(0.5853981633974485,pi/2) q37[9];
cx q37[9],q37[8];
rz(0.5853981633974479) q37[8];
sx q37[9];
cx q37[9],q37[8];
rz(-pi) q37[8];
u3(pi,-3*pi/4,3*pi/4) q37[9];
barrier q37[0],q37[1],q37[2],q37[3],q37[4],q37[5],q37[6],q37[7],q37[8],q37[9],q37[10],q37[11],q37[12],q37[13],q37[14],q37[15],q37[16],q37[17],q37[18],q37[19],q37[20],q37[21],q37[22],q37[23];
u2(-pi,-pi/2) q37[0];
sx q37[1];
sx q37[10];
u2(-pi,-pi/2) q37[11];
cx q37[11],q37[10];
rz(-0.9853981633974493) q37[10];
u2(0.9853981633974493,pi/2) q37[11];
cx q37[11],q37[10];
rz(0.9853981633974496) q37[10];
sx q37[11];
cx q37[11],q37[10];
u2(-pi/2,-pi/2) q37[10];
u2(pi/2,0) q37[11];
cx q37[0],q37[11];
u2(0.8000000000000007,pi/2) q37[0];
rz(-0.8000000000000006) q37[11];
cx q37[0],q37[11];
sx q37[0];
rz(0.8000000000000008) q37[11];
cx q37[0],q37[11];
u3(pi,-3*pi/4,3*pi/4) q37[0];
u2(0,pi/2) q37[11];
cx q37[11],q37[10];
rz(0.5853981633974487) q37[10];
u2(0.5853981633974485,pi/2) q37[11];
cx q37[11],q37[10];
rz(0.5853981633974479) q37[10];
sx q37[11];
cx q37[11],q37[10];
rz(-pi) q37[10];
u3(pi,-3*pi/4,3*pi/4) q37[11];
u2(-pi,-pi/2) q37[12];
cx q37[12],q37[1];
rz(-0.9853981633974493) q37[1];
u2(0.9853981633974493,pi/2) q37[12];
cx q37[12],q37[1];
rz(0.9853981633974496) q37[1];
sx q37[12];
cx q37[12],q37[1];
u2(-pi/2,-pi/2) q37[1];
u2(pi/2,0) q37[12];
u2(-pi,-pi/2) q37[13];
cx q37[13],q37[12];
rz(-0.8000000000000006) q37[12];
u2(0.8000000000000007,pi/2) q37[13];
cx q37[13],q37[12];
rz(0.8000000000000008) q37[12];
sx q37[13];
cx q37[13],q37[12];
u2(0,pi/2) q37[12];
cx q37[12],q37[1];
rz(0.5853981633974487) q37[1];
u2(0.5853981633974485,pi/2) q37[12];
cx q37[12],q37[1];
rz(0.5853981633974479) q37[1];
sx q37[12];
cx q37[12],q37[1];
rz(-pi) q37[1];
u3(pi,-3*pi/4,3*pi/4) q37[12];
u3(pi,-3*pi/4,3*pi/4) q37[13];
u2(-pi,-pi/2) q37[16];
u2(-pi,-pi/2) q37[17];
sx q37[2];
u2(-pi,-pi/2) q37[20];
u2(-pi,-pi/2) q37[21];
u2(-pi,-pi/2) q37[3];
cx q37[3],q37[2];
rz(-0.9853981633974493) q37[2];
u2(0.9853981633974493,pi/2) q37[3];
cx q37[3],q37[2];
rz(0.9853981633974496) q37[2];
sx q37[3];
cx q37[3],q37[2];
u2(-pi/2,-pi/2) q37[2];
u2(pi/2,0) q37[3];
u2(-pi,-pi/2) q37[4];
cx q37[4],q37[3];
rz(-0.8000000000000006) q37[3];
u2(0.8000000000000007,pi/2) q37[4];
cx q37[4],q37[3];
rz(0.8000000000000008) q37[3];
sx q37[4];
cx q37[4],q37[3];
u2(0,pi/2) q37[3];
cx q37[3],q37[2];
rz(0.5853981633974487) q37[2];
u2(0.5853981633974485,pi/2) q37[3];
cx q37[3],q37[2];
rz(0.5853981633974479) q37[2];
sx q37[3];
cx q37[3],q37[2];
rz(-pi) q37[2];
u3(pi,-3*pi/4,3*pi/4) q37[3];
u3(pi,-3*pi/4,3*pi/4) q37[4];
sx q37[5];
cx q37[16],q37[5];
u2(0.9853981633974493,pi/2) q37[16];
rz(-0.9853981633974493) q37[5];
cx q37[16],q37[5];
sx q37[16];
rz(0.9853981633974496) q37[5];
cx q37[16],q37[5];
u2(pi/2,0) q37[16];
cx q37[17],q37[16];
rz(-0.8000000000000006) q37[16];
u2(0.8000000000000007,pi/2) q37[17];
cx q37[17],q37[16];
rz(0.8000000000000008) q37[16];
sx q37[17];
cx q37[17],q37[16];
u2(0,pi/2) q37[16];
u3(pi,-3*pi/4,3*pi/4) q37[17];
u2(-pi/2,-pi/2) q37[5];
cx q37[16],q37[5];
u2(0.5853981633974485,pi/2) q37[16];
rz(0.5853981633974487) q37[5];
cx q37[16],q37[5];
sx q37[16];
rz(0.5853981633974479) q37[5];
cx q37[16],q37[5];
u3(pi,-3*pi/4,3*pi/4) q37[16];
rz(-pi) q37[5];
sx q37[6];
u2(-pi,-pi/2) q37[7];
cx q37[7],q37[6];
rz(-0.9853981633974493) q37[6];
u2(0.9853981633974493,pi/2) q37[7];
cx q37[7],q37[6];
rz(0.9853981633974496) q37[6];
sx q37[7];
cx q37[7],q37[6];
u2(-pi/2,-pi/2) q37[6];
u2(pi/2,0) q37[7];
u2(-pi,-pi/2) q37[8];
cx q37[8],q37[7];
rz(-0.8000000000000006) q37[7];
u2(0.8000000000000007,pi/2) q37[8];
cx q37[8],q37[7];
rz(0.8000000000000008) q37[7];
sx q37[8];
cx q37[8],q37[7];
u2(0,pi/2) q37[7];
cx q37[7],q37[6];
rz(0.5853981633974487) q37[6];
u2(0.5853981633974485,pi/2) q37[7];
cx q37[7],q37[6];
rz(0.5853981633974479) q37[6];
sx q37[7];
cx q37[7],q37[6];
rz(-pi) q37[6];
u3(pi,-3*pi/4,3*pi/4) q37[7];
u3(pi,-3*pi/4,3*pi/4) q37[8];
sx q37[9];
cx q37[20],q37[9];
u2(0.9853981633974493,pi/2) q37[20];
rz(-0.9853981633974493) q37[9];
cx q37[20],q37[9];
sx q37[20];
rz(0.9853981633974496) q37[9];
cx q37[20],q37[9];
u2(pi/2,0) q37[20];
cx q37[21],q37[20];
rz(-0.8000000000000006) q37[20];
u2(0.8000000000000007,pi/2) q37[21];
cx q37[21],q37[20];
rz(0.8000000000000008) q37[20];
sx q37[21];
cx q37[21],q37[20];
u2(0,pi/2) q37[20];
u3(pi,-3*pi/4,3*pi/4) q37[21];
u2(-pi/2,-pi/2) q37[9];
cx q37[20],q37[9];
u2(0.5853981633974485,pi/2) q37[20];
rz(0.5853981633974487) q37[9];
cx q37[20],q37[9];
sx q37[20];
rz(0.5853981633974479) q37[9];
cx q37[20],q37[9];
u3(pi,-3*pi/4,3*pi/4) q37[20];
rz(-pi) q37[9];
barrier q37[0],q37[1],q37[2],q37[3],q37[4],q37[5],q37[6],q37[7],q37[8],q37[9],q37[10],q37[11],q37[12],q37[13],q37[14],q37[15],q37[16],q37[17],q37[18],q37[19],q37[20],q37[21],q37[22],q37[23];
sx q37[0];
u2(-pi,-pi/2) q37[1];
cx q37[1],q37[0];
rz(-0.9853981633974493) q37[0];
u2(0.9853981633974493,pi/2) q37[1];
cx q37[1],q37[0];
rz(0.9853981633974496) q37[0];
sx q37[1];
cx q37[1],q37[0];
u2(-pi/2,-pi/2) q37[0];
u2(pi/2,0) q37[1];
u2(-pi,-pi/2) q37[10];
sx q37[11];
u2(-pi,-pi/2) q37[14];
u2(-pi,-pi/2) q37[15];
u2(-pi,-pi/2) q37[18];
u2(-pi,-pi/2) q37[19];
u2(-pi,-pi/2) q37[2];
cx q37[2],q37[1];
rz(-0.8000000000000006) q37[1];
u2(0.8000000000000007,pi/2) q37[2];
cx q37[2],q37[1];
rz(0.8000000000000008) q37[1];
sx q37[2];
cx q37[2],q37[1];
u2(0,pi/2) q37[1];
cx q37[1],q37[0];
rz(0.5853981633974487) q37[0];
u2(0.5853981633974485,pi/2) q37[1];
cx q37[1],q37[0];
rz(0.5853981633974479) q37[0];
sx q37[1];
cx q37[1],q37[0];
rz(-pi) q37[0];
u3(pi,-3*pi/4,3*pi/4) q37[1];
u3(pi,-3*pi/4,3*pi/4) q37[2];
u2(-pi,-pi/2) q37[22];
cx q37[22],q37[11];
rz(-0.9853981633974493) q37[11];
u2(0.9853981633974493,pi/2) q37[22];
cx q37[22],q37[11];
rz(0.9853981633974496) q37[11];
sx q37[22];
cx q37[22],q37[11];
u2(-pi/2,-pi/2) q37[11];
u2(pi/2,0) q37[22];
u2(-pi,-pi/2) q37[23];
cx q37[23],q37[22];
rz(-0.8000000000000006) q37[22];
u2(0.8000000000000007,pi/2) q37[23];
cx q37[23],q37[22];
rz(0.8000000000000008) q37[22];
sx q37[23];
cx q37[23],q37[22];
u2(0,pi/2) q37[22];
cx q37[22],q37[11];
rz(0.5853981633974487) q37[11];
u2(0.5853981633974485,pi/2) q37[22];
cx q37[22],q37[11];
rz(0.5853981633974479) q37[11];
sx q37[22];
cx q37[22],q37[11];
rz(-pi) q37[11];
u3(pi,-3*pi/4,3*pi/4) q37[22];
u3(pi,-3*pi/4,3*pi/4) q37[23];
sx q37[3];
cx q37[14],q37[3];
u2(0.9853981633974493,pi/2) q37[14];
rz(-0.9853981633974493) q37[3];
cx q37[14],q37[3];
sx q37[14];
rz(0.9853981633974496) q37[3];
cx q37[14],q37[3];
u2(pi/2,0) q37[14];
cx q37[15],q37[14];
rz(-0.8000000000000006) q37[14];
u2(0.8000000000000007,pi/2) q37[15];
cx q37[15],q37[14];
rz(0.8000000000000008) q37[14];
sx q37[15];
cx q37[15],q37[14];
u2(0,pi/2) q37[14];
u3(pi,-3*pi/4,3*pi/4) q37[15];
u2(-pi/2,-pi/2) q37[3];
cx q37[14],q37[3];
u2(0.5853981633974485,pi/2) q37[14];
rz(0.5853981633974487) q37[3];
cx q37[14],q37[3];
sx q37[14];
rz(0.5853981633974479) q37[3];
cx q37[14],q37[3];
u3(pi,-3*pi/4,3*pi/4) q37[14];
rz(-pi) q37[3];
sx q37[4];
u2(-pi,-pi/2) q37[5];
cx q37[5],q37[4];
rz(-0.9853981633974493) q37[4];
u2(0.9853981633974493,pi/2) q37[5];
cx q37[5],q37[4];
rz(0.9853981633974496) q37[4];
sx q37[5];
cx q37[5],q37[4];
u2(-pi/2,-pi/2) q37[4];
u2(pi/2,0) q37[5];
u2(-pi,-pi/2) q37[6];
cx q37[6],q37[5];
rz(-0.8000000000000006) q37[5];
u2(0.8000000000000007,pi/2) q37[6];
cx q37[6],q37[5];
rz(0.8000000000000008) q37[5];
sx q37[6];
cx q37[6],q37[5];
u2(0,pi/2) q37[5];
cx q37[5],q37[4];
rz(0.5853981633974487) q37[4];
u2(0.5853981633974485,pi/2) q37[5];
cx q37[5],q37[4];
rz(0.5853981633974479) q37[4];
sx q37[5];
cx q37[5],q37[4];
rz(-pi) q37[4];
u3(pi,-3*pi/4,3*pi/4) q37[5];
u3(pi,-3*pi/4,3*pi/4) q37[6];
sx q37[7];
cx q37[18],q37[7];
u2(0.9853981633974493,pi/2) q37[18];
rz(-0.9853981633974493) q37[7];
cx q37[18],q37[7];
sx q37[18];
rz(0.9853981633974496) q37[7];
cx q37[18],q37[7];
u2(pi/2,0) q37[18];
cx q37[19],q37[18];
rz(-0.8000000000000006) q37[18];
u2(0.8000000000000007,pi/2) q37[19];
cx q37[19],q37[18];
rz(0.8000000000000008) q37[18];
sx q37[19];
cx q37[19],q37[18];
u2(0,pi/2) q37[18];
u3(pi,-3*pi/4,3*pi/4) q37[19];
u2(-pi/2,-pi/2) q37[7];
cx q37[18],q37[7];
u2(0.5853981633974485,pi/2) q37[18];
rz(0.5853981633974487) q37[7];
cx q37[18],q37[7];
sx q37[18];
rz(0.5853981633974479) q37[7];
cx q37[18],q37[7];
u3(pi,-3*pi/4,3*pi/4) q37[18];
rz(-pi) q37[7];
sx q37[8];
u2(-pi,-pi/2) q37[9];
cx q37[9],q37[8];
rz(-0.9853981633974493) q37[8];
u2(0.9853981633974493,pi/2) q37[9];
cx q37[9],q37[8];
rz(0.9853981633974496) q37[8];
sx q37[9];
cx q37[9],q37[8];
u2(-pi/2,-pi/2) q37[8];
u2(pi/2,0) q37[9];
cx q37[10],q37[9];
u2(0.8000000000000007,pi/2) q37[10];
rz(-0.8000000000000006) q37[9];
cx q37[10],q37[9];
sx q37[10];
rz(0.8000000000000008) q37[9];
cx q37[10],q37[9];
u3(pi,-3*pi/4,3*pi/4) q37[10];
u2(0,pi/2) q37[9];
cx q37[9],q37[8];
rz(0.5853981633974487) q37[8];
u2(0.5853981633974485,pi/2) q37[9];
cx q37[9],q37[8];
rz(0.5853981633974479) q37[8];
sx q37[9];
cx q37[9],q37[8];
rz(-pi) q37[8];
u3(pi,-3*pi/4,3*pi/4) q37[9];
barrier q37[0],q37[1],q37[2],q37[3],q37[4],q37[5],q37[6],q37[7],q37[8],q37[9],q37[10],q37[11],q37[12],q37[13],q37[14],q37[15],q37[16],q37[17],q37[18],q37[19],q37[20],q37[21],q37[22],q37[23];
u2(-pi,-pi/2) q37[0];
sx q37[1];
sx q37[10];
u2(-pi,-pi/2) q37[11];
cx q37[11],q37[10];
rz(-0.9853981633974493) q37[10];
u2(0.9853981633974493,pi/2) q37[11];
cx q37[11],q37[10];
rz(0.9853981633974496) q37[10];
sx q37[11];
cx q37[11],q37[10];
u2(-pi/2,-pi/2) q37[10];
u2(pi/2,0) q37[11];
cx q37[0],q37[11];
u2(0.8000000000000007,pi/2) q37[0];
rz(-0.8000000000000006) q37[11];
cx q37[0],q37[11];
sx q37[0];
rz(0.8000000000000008) q37[11];
cx q37[0],q37[11];
u3(pi,-3*pi/4,3*pi/4) q37[0];
u2(0,pi/2) q37[11];
cx q37[11],q37[10];
rz(0.5853981633974487) q37[10];
u2(0.5853981633974485,pi/2) q37[11];
cx q37[11],q37[10];
rz(0.5853981633974479) q37[10];
sx q37[11];
cx q37[11],q37[10];
rz(-pi) q37[10];
u3(pi,-3*pi/4,3*pi/4) q37[11];
u2(-pi,-pi/2) q37[12];
cx q37[12],q37[1];
rz(-0.9853981633974493) q37[1];
u2(0.9853981633974493,pi/2) q37[12];
cx q37[12],q37[1];
rz(0.9853981633974496) q37[1];
sx q37[12];
cx q37[12],q37[1];
u2(-pi/2,-pi/2) q37[1];
u2(pi/2,0) q37[12];
u2(-pi,-pi/2) q37[13];
cx q37[13],q37[12];
rz(-0.8000000000000006) q37[12];
u2(0.8000000000000007,pi/2) q37[13];
cx q37[13],q37[12];
rz(0.8000000000000008) q37[12];
sx q37[13];
cx q37[13],q37[12];
u2(0,pi/2) q37[12];
cx q37[12],q37[1];
rz(0.5853981633974487) q37[1];
u2(0.5853981633974485,pi/2) q37[12];
cx q37[12],q37[1];
rz(0.5853981633974479) q37[1];
sx q37[12];
cx q37[12],q37[1];
rz(-pi) q37[1];
u3(pi,-3*pi/4,3*pi/4) q37[12];
u3(pi,-3*pi/4,3*pi/4) q37[13];
u2(-pi,-pi/2) q37[16];
u2(-pi,-pi/2) q37[17];
sx q37[2];
u2(-pi,-pi/2) q37[20];
u2(-pi,-pi/2) q37[21];
u2(-pi,-pi/2) q37[3];
cx q37[3],q37[2];
rz(-0.9853981633974493) q37[2];
u2(0.9853981633974493,pi/2) q37[3];
cx q37[3],q37[2];
rz(0.9853981633974496) q37[2];
sx q37[3];
cx q37[3],q37[2];
u2(-pi/2,-pi/2) q37[2];
u2(pi/2,0) q37[3];
u2(-pi,-pi/2) q37[4];
cx q37[4],q37[3];
rz(-0.8000000000000006) q37[3];
u2(0.8000000000000007,pi/2) q37[4];
cx q37[4],q37[3];
rz(0.8000000000000008) q37[3];
sx q37[4];
cx q37[4],q37[3];
u2(0,pi/2) q37[3];
cx q37[3],q37[2];
rz(0.5853981633974487) q37[2];
u2(0.5853981633974485,pi/2) q37[3];
cx q37[3],q37[2];
rz(0.5853981633974479) q37[2];
sx q37[3];
cx q37[3],q37[2];
rz(-pi) q37[2];
u3(pi,-3*pi/4,3*pi/4) q37[3];
u3(pi,-3*pi/4,3*pi/4) q37[4];
sx q37[5];
cx q37[16],q37[5];
u2(0.9853981633974493,pi/2) q37[16];
rz(-0.9853981633974493) q37[5];
cx q37[16],q37[5];
sx q37[16];
rz(0.9853981633974496) q37[5];
cx q37[16],q37[5];
u2(pi/2,0) q37[16];
cx q37[17],q37[16];
rz(-0.8000000000000006) q37[16];
u2(0.8000000000000007,pi/2) q37[17];
cx q37[17],q37[16];
rz(0.8000000000000008) q37[16];
sx q37[17];
cx q37[17],q37[16];
u2(0,pi/2) q37[16];
u3(pi,-3*pi/4,3*pi/4) q37[17];
u2(-pi/2,-pi/2) q37[5];
cx q37[16],q37[5];
u2(0.5853981633974485,pi/2) q37[16];
rz(0.5853981633974487) q37[5];
cx q37[16],q37[5];
sx q37[16];
rz(0.5853981633974479) q37[5];
cx q37[16],q37[5];
u3(pi,-3*pi/4,3*pi/4) q37[16];
rz(-pi) q37[5];
sx q37[6];
u2(-pi,-pi/2) q37[7];
cx q37[7],q37[6];
rz(-0.9853981633974493) q37[6];
u2(0.9853981633974493,pi/2) q37[7];
cx q37[7],q37[6];
rz(0.9853981633974496) q37[6];
sx q37[7];
cx q37[7],q37[6];
u2(-pi/2,-pi/2) q37[6];
u2(pi/2,0) q37[7];
u2(-pi,-pi/2) q37[8];
cx q37[8],q37[7];
rz(-0.8000000000000006) q37[7];
u2(0.8000000000000007,pi/2) q37[8];
cx q37[8],q37[7];
rz(0.8000000000000008) q37[7];
sx q37[8];
cx q37[8],q37[7];
u2(0,pi/2) q37[7];
cx q37[7],q37[6];
rz(0.5853981633974487) q37[6];
u2(0.5853981633974485,pi/2) q37[7];
cx q37[7],q37[6];
rz(0.5853981633974479) q37[6];
sx q37[7];
cx q37[7],q37[6];
rz(-pi) q37[6];
u3(pi,-3*pi/4,3*pi/4) q37[7];
u3(pi,-3*pi/4,3*pi/4) q37[8];
sx q37[9];
cx q37[20],q37[9];
u2(0.9853981633974493,pi/2) q37[20];
rz(-0.9853981633974493) q37[9];
cx q37[20],q37[9];
sx q37[20];
rz(0.9853981633974496) q37[9];
cx q37[20],q37[9];
u2(pi/2,0) q37[20];
cx q37[21],q37[20];
rz(-0.8000000000000006) q37[20];
u2(0.8000000000000007,pi/2) q37[21];
cx q37[21],q37[20];
rz(0.8000000000000008) q37[20];
sx q37[21];
cx q37[21],q37[20];
u2(0,pi/2) q37[20];
u3(pi,-3*pi/4,3*pi/4) q37[21];
u2(-pi/2,-pi/2) q37[9];
cx q37[20],q37[9];
u2(0.5853981633974485,pi/2) q37[20];
rz(0.5853981633974487) q37[9];
cx q37[20],q37[9];
sx q37[20];
rz(0.5853981633974479) q37[9];
cx q37[20],q37[9];
u3(pi,-3*pi/4,3*pi/4) q37[20];
rz(-pi) q37[9];
barrier q37[0],q37[1],q37[2],q37[3],q37[4],q37[5],q37[6],q37[7],q37[8],q37[9],q37[10],q37[11],q37[12],q37[13],q37[14],q37[15],q37[16],q37[17],q37[18],q37[19],q37[20],q37[21],q37[22],q37[23];
sx q37[0];
u2(-pi,-pi/2) q37[1];
cx q37[1],q37[0];
rz(-0.9853981633974493) q37[0];
u2(0.9853981633974493,pi/2) q37[1];
cx q37[1],q37[0];
rz(0.9853981633974496) q37[0];
sx q37[1];
cx q37[1],q37[0];
u2(-pi/2,-pi/2) q37[0];
u2(pi/2,0) q37[1];
u2(-pi,-pi/2) q37[10];
sx q37[11];
u2(-pi,-pi/2) q37[14];
u2(-pi,-pi/2) q37[15];
u2(-pi,-pi/2) q37[18];
u2(-pi,-pi/2) q37[19];
u2(-pi,-pi/2) q37[2];
cx q37[2],q37[1];
rz(-0.8000000000000006) q37[1];
u2(0.8000000000000007,pi/2) q37[2];
cx q37[2],q37[1];
rz(0.8000000000000008) q37[1];
sx q37[2];
cx q37[2],q37[1];
u2(0,pi/2) q37[1];
cx q37[1],q37[0];
rz(0.5853981633974487) q37[0];
u2(0.5853981633974485,pi/2) q37[1];
cx q37[1],q37[0];
rz(0.5853981633974479) q37[0];
sx q37[1];
cx q37[1],q37[0];
rz(-pi) q37[0];
u3(pi,-3*pi/4,3*pi/4) q37[1];
u3(pi,-3*pi/4,3*pi/4) q37[2];
u2(-pi,-pi/2) q37[22];
cx q37[22],q37[11];
rz(-0.9853981633974493) q37[11];
u2(0.9853981633974493,pi/2) q37[22];
cx q37[22],q37[11];
rz(0.9853981633974496) q37[11];
sx q37[22];
cx q37[22],q37[11];
u2(-pi/2,-pi/2) q37[11];
u2(pi/2,0) q37[22];
u2(-pi,-pi/2) q37[23];
cx q37[23],q37[22];
rz(-0.8000000000000006) q37[22];
u2(0.8000000000000007,pi/2) q37[23];
cx q37[23],q37[22];
rz(0.8000000000000008) q37[22];
sx q37[23];
cx q37[23],q37[22];
u2(0,pi/2) q37[22];
cx q37[22],q37[11];
rz(0.5853981633974487) q37[11];
u2(0.5853981633974485,pi/2) q37[22];
cx q37[22],q37[11];
rz(0.5853981633974479) q37[11];
sx q37[22];
cx q37[22],q37[11];
rz(-pi) q37[11];
u3(pi,-3*pi/4,3*pi/4) q37[22];
u3(pi,-3*pi/4,3*pi/4) q37[23];
sx q37[3];
cx q37[14],q37[3];
u2(0.9853981633974493,pi/2) q37[14];
rz(-0.9853981633974493) q37[3];
cx q37[14],q37[3];
sx q37[14];
rz(0.9853981633974496) q37[3];
cx q37[14],q37[3];
u2(pi/2,0) q37[14];
cx q37[15],q37[14];
rz(-0.8000000000000006) q37[14];
u2(0.8000000000000007,pi/2) q37[15];
cx q37[15],q37[14];
rz(0.8000000000000008) q37[14];
sx q37[15];
cx q37[15],q37[14];
u2(0,pi/2) q37[14];
u3(pi,-3*pi/4,3*pi/4) q37[15];
u2(-pi/2,-pi/2) q37[3];
cx q37[14],q37[3];
u2(0.5853981633974485,pi/2) q37[14];
rz(0.5853981633974487) q37[3];
cx q37[14],q37[3];
sx q37[14];
rz(0.5853981633974479) q37[3];
cx q37[14],q37[3];
u3(pi,-3*pi/4,3*pi/4) q37[14];
rz(-pi) q37[3];
sx q37[4];
u2(-pi,-pi/2) q37[5];
cx q37[5],q37[4];
rz(-0.9853981633974493) q37[4];
u2(0.9853981633974493,pi/2) q37[5];
cx q37[5],q37[4];
rz(0.9853981633974496) q37[4];
sx q37[5];
cx q37[5],q37[4];
u2(-pi/2,-pi/2) q37[4];
u2(pi/2,0) q37[5];
u2(-pi,-pi/2) q37[6];
cx q37[6],q37[5];
rz(-0.8000000000000006) q37[5];
u2(0.8000000000000007,pi/2) q37[6];
cx q37[6],q37[5];
rz(0.8000000000000008) q37[5];
sx q37[6];
cx q37[6],q37[5];
u2(0,pi/2) q37[5];
cx q37[5],q37[4];
rz(0.5853981633974487) q37[4];
u2(0.5853981633974485,pi/2) q37[5];
cx q37[5],q37[4];
rz(0.5853981633974479) q37[4];
sx q37[5];
cx q37[5],q37[4];
rz(-pi) q37[4];
u3(pi,-3*pi/4,3*pi/4) q37[5];
u3(pi,-3*pi/4,3*pi/4) q37[6];
sx q37[7];
cx q37[18],q37[7];
u2(0.9853981633974493,pi/2) q37[18];
rz(-0.9853981633974493) q37[7];
cx q37[18],q37[7];
sx q37[18];
rz(0.9853981633974496) q37[7];
cx q37[18],q37[7];
u2(pi/2,0) q37[18];
cx q37[19],q37[18];
rz(-0.8000000000000006) q37[18];
u2(0.8000000000000007,pi/2) q37[19];
cx q37[19],q37[18];
rz(0.8000000000000008) q37[18];
sx q37[19];
cx q37[19],q37[18];
u2(0,pi/2) q37[18];
u3(pi,-3*pi/4,3*pi/4) q37[19];
u2(-pi/2,-pi/2) q37[7];
cx q37[18],q37[7];
u2(0.5853981633974485,pi/2) q37[18];
rz(0.5853981633974487) q37[7];
cx q37[18],q37[7];
sx q37[18];
rz(0.5853981633974479) q37[7];
cx q37[18],q37[7];
u3(pi,-3*pi/4,3*pi/4) q37[18];
rz(-pi) q37[7];
sx q37[8];
u2(-pi,-pi/2) q37[9];
cx q37[9],q37[8];
rz(-0.9853981633974493) q37[8];
u2(0.9853981633974493,pi/2) q37[9];
cx q37[9],q37[8];
rz(0.9853981633974496) q37[8];
sx q37[9];
cx q37[9],q37[8];
u2(-pi/2,-pi/2) q37[8];
u2(pi/2,0) q37[9];
cx q37[10],q37[9];
u2(0.8000000000000007,pi/2) q37[10];
rz(-0.8000000000000006) q37[9];
cx q37[10],q37[9];
sx q37[10];
rz(0.8000000000000008) q37[9];
cx q37[10],q37[9];
u3(pi,-3*pi/4,3*pi/4) q37[10];
u2(0,pi/2) q37[9];
cx q37[9],q37[8];
rz(0.5853981633974487) q37[8];
u2(0.5853981633974485,pi/2) q37[9];
cx q37[9],q37[8];
rz(0.5853981633974479) q37[8];
sx q37[9];
cx q37[9],q37[8];
rz(-pi) q37[8];
u3(pi,-3*pi/4,3*pi/4) q37[9];
barrier q37[0],q37[1],q37[2],q37[3],q37[4],q37[5],q37[6],q37[7],q37[8],q37[9],q37[10],q37[11],q37[12],q37[13],q37[14],q37[15],q37[16],q37[17],q37[18],q37[19],q37[20],q37[21],q37[22],q37[23];
u2(-pi,-pi/2) q37[0];
sx q37[1];
sx q37[10];
u2(-pi,-pi/2) q37[11];
cx q37[11],q37[10];
rz(-0.9853981633974493) q37[10];
u2(0.9853981633974493,pi/2) q37[11];
cx q37[11],q37[10];
rz(0.9853981633974496) q37[10];
sx q37[11];
cx q37[11],q37[10];
u2(-pi/2,-pi/2) q37[10];
u2(pi/2,0) q37[11];
cx q37[0],q37[11];
u2(0.8000000000000007,pi/2) q37[0];
rz(-0.8000000000000006) q37[11];
cx q37[0],q37[11];
sx q37[0];
rz(0.8000000000000008) q37[11];
cx q37[0],q37[11];
u3(pi,-3*pi/4,3*pi/4) q37[0];
u2(0,pi/2) q37[11];
cx q37[11],q37[10];
rz(0.5853981633974487) q37[10];
u2(0.5853981633974485,pi/2) q37[11];
cx q37[11],q37[10];
rz(0.5853981633974479) q37[10];
sx q37[11];
cx q37[11],q37[10];
rz(-pi) q37[10];
u3(pi,-3*pi/4,3*pi/4) q37[11];
u2(-pi,-pi/2) q37[12];
cx q37[12],q37[1];
rz(-0.9853981633974493) q37[1];
u2(0.9853981633974493,pi/2) q37[12];
cx q37[12],q37[1];
rz(0.9853981633974496) q37[1];
sx q37[12];
cx q37[12],q37[1];
u2(-pi/2,-pi/2) q37[1];
u2(pi/2,0) q37[12];
u2(-pi,-pi/2) q37[13];
cx q37[13],q37[12];
rz(-0.8000000000000006) q37[12];
u2(0.8000000000000007,pi/2) q37[13];
cx q37[13],q37[12];
rz(0.8000000000000008) q37[12];
sx q37[13];
cx q37[13],q37[12];
u2(0,pi/2) q37[12];
cx q37[12],q37[1];
rz(0.5853981633974487) q37[1];
u2(0.5853981633974485,pi/2) q37[12];
cx q37[12],q37[1];
rz(0.5853981633974479) q37[1];
sx q37[12];
cx q37[12],q37[1];
rz(-pi) q37[1];
u3(pi,-3*pi/4,3*pi/4) q37[12];
u3(pi,-3*pi/4,3*pi/4) q37[13];
u2(-pi,-pi/2) q37[16];
u2(-pi,-pi/2) q37[17];
sx q37[2];
u2(-pi,-pi/2) q37[20];
u2(-pi,-pi/2) q37[21];
u2(-pi,-pi/2) q37[3];
cx q37[3],q37[2];
rz(-0.9853981633974493) q37[2];
u2(0.9853981633974493,pi/2) q37[3];
cx q37[3],q37[2];
rz(0.9853981633974496) q37[2];
sx q37[3];
cx q37[3],q37[2];
u2(-pi/2,-pi/2) q37[2];
u2(pi/2,0) q37[3];
u2(-pi,-pi/2) q37[4];
cx q37[4],q37[3];
rz(-0.8000000000000006) q37[3];
u2(0.8000000000000007,pi/2) q37[4];
cx q37[4],q37[3];
rz(0.8000000000000008) q37[3];
sx q37[4];
cx q37[4],q37[3];
u2(0,pi/2) q37[3];
cx q37[3],q37[2];
rz(0.5853981633974487) q37[2];
u2(0.5853981633974485,pi/2) q37[3];
cx q37[3],q37[2];
rz(0.5853981633974479) q37[2];
sx q37[3];
cx q37[3],q37[2];
rz(-pi) q37[2];
u3(pi,-3*pi/4,3*pi/4) q37[3];
u3(pi,-3*pi/4,3*pi/4) q37[4];
sx q37[5];
cx q37[16],q37[5];
u2(0.9853981633974493,pi/2) q37[16];
rz(-0.9853981633974493) q37[5];
cx q37[16],q37[5];
sx q37[16];
rz(0.9853981633974496) q37[5];
cx q37[16],q37[5];
u2(pi/2,0) q37[16];
cx q37[17],q37[16];
rz(-0.8000000000000006) q37[16];
u2(0.8000000000000007,pi/2) q37[17];
cx q37[17],q37[16];
rz(0.8000000000000008) q37[16];
sx q37[17];
cx q37[17],q37[16];
u2(0,pi/2) q37[16];
u3(pi,-3*pi/4,3*pi/4) q37[17];
u2(-pi/2,-pi/2) q37[5];
cx q37[16],q37[5];
u2(0.5853981633974485,pi/2) q37[16];
rz(0.5853981633974487) q37[5];
cx q37[16],q37[5];
sx q37[16];
rz(0.5853981633974479) q37[5];
cx q37[16],q37[5];
u3(pi,-3*pi/4,3*pi/4) q37[16];
rz(-pi) q37[5];
sx q37[6];
u2(-pi,-pi/2) q37[7];
cx q37[7],q37[6];
rz(-0.9853981633974493) q37[6];
u2(0.9853981633974493,pi/2) q37[7];
cx q37[7],q37[6];
rz(0.9853981633974496) q37[6];
sx q37[7];
cx q37[7],q37[6];
u2(-pi/2,-pi/2) q37[6];
u2(pi/2,0) q37[7];
u2(-pi,-pi/2) q37[8];
cx q37[8],q37[7];
rz(-0.8000000000000006) q37[7];
u2(0.8000000000000007,pi/2) q37[8];
cx q37[8],q37[7];
rz(0.8000000000000008) q37[7];
sx q37[8];
cx q37[8],q37[7];
u2(0,pi/2) q37[7];
cx q37[7],q37[6];
rz(0.5853981633974487) q37[6];
u2(0.5853981633974485,pi/2) q37[7];
cx q37[7],q37[6];
rz(0.5853981633974479) q37[6];
sx q37[7];
cx q37[7],q37[6];
rz(-pi) q37[6];
u3(pi,-3*pi/4,3*pi/4) q37[7];
u3(pi,-3*pi/4,3*pi/4) q37[8];
sx q37[9];
cx q37[20],q37[9];
u2(0.9853981633974493,pi/2) q37[20];
rz(-0.9853981633974493) q37[9];
cx q37[20],q37[9];
sx q37[20];
rz(0.9853981633974496) q37[9];
cx q37[20],q37[9];
u2(pi/2,0) q37[20];
cx q37[21],q37[20];
rz(-0.8000000000000006) q37[20];
u2(0.8000000000000007,pi/2) q37[21];
cx q37[21],q37[20];
rz(0.8000000000000008) q37[20];
sx q37[21];
cx q37[21],q37[20];
u2(0,pi/2) q37[20];
u3(pi,-3*pi/4,3*pi/4) q37[21];
u2(-pi/2,-pi/2) q37[9];
cx q37[20],q37[9];
u2(0.5853981633974485,pi/2) q37[20];
rz(0.5853981633974487) q37[9];
cx q37[20],q37[9];
sx q37[20];
rz(0.5853981633974479) q37[9];
cx q37[20],q37[9];
u3(pi,-3*pi/4,3*pi/4) q37[20];
rz(-pi) q37[9];
barrier q37[0],q37[1],q37[2],q37[3],q37[4],q37[5],q37[6],q37[7],q37[8],q37[9],q37[10],q37[11],q37[12],q37[13],q37[14],q37[15],q37[16],q37[17],q37[18],q37[19],q37[20],q37[21],q37[22],q37[23];
sx q37[0];
u2(-pi,-pi/2) q37[1];
cx q37[1],q37[0];
rz(-0.9853981633974493) q37[0];
u2(0.9853981633974493,pi/2) q37[1];
cx q37[1],q37[0];
rz(0.9853981633974496) q37[0];
sx q37[1];
cx q37[1],q37[0];
u2(-pi/2,-pi/2) q37[0];
u2(pi/2,0) q37[1];
u2(-pi,-pi/2) q37[10];
sx q37[11];
u2(-pi,-pi/2) q37[14];
u2(-pi,-pi/2) q37[15];
u2(-pi,-pi/2) q37[18];
u2(-pi,-pi/2) q37[19];
u2(-pi,-pi/2) q37[2];
cx q37[2],q37[1];
rz(-0.8000000000000006) q37[1];
u2(0.8000000000000007,pi/2) q37[2];
cx q37[2],q37[1];
rz(0.8000000000000008) q37[1];
sx q37[2];
cx q37[2],q37[1];
u2(0,pi/2) q37[1];
cx q37[1],q37[0];
rz(0.5853981633974487) q37[0];
u2(0.5853981633974485,pi/2) q37[1];
cx q37[1],q37[0];
rz(0.5853981633974479) q37[0];
sx q37[1];
cx q37[1],q37[0];
rz(-pi) q37[0];
u3(pi,-3*pi/4,3*pi/4) q37[1];
u3(pi,-3*pi/4,3*pi/4) q37[2];
u2(-pi,-pi/2) q37[22];
cx q37[22],q37[11];
rz(-0.9853981633974493) q37[11];
u2(0.9853981633974493,pi/2) q37[22];
cx q37[22],q37[11];
rz(0.9853981633974496) q37[11];
sx q37[22];
cx q37[22],q37[11];
u2(-pi/2,-pi/2) q37[11];
u2(pi/2,0) q37[22];
u2(-pi,-pi/2) q37[23];
cx q37[23],q37[22];
rz(-0.8000000000000006) q37[22];
u2(0.8000000000000007,pi/2) q37[23];
cx q37[23],q37[22];
rz(0.8000000000000008) q37[22];
sx q37[23];
cx q37[23],q37[22];
u2(0,pi/2) q37[22];
cx q37[22],q37[11];
rz(0.5853981633974487) q37[11];
u2(0.5853981633974485,pi/2) q37[22];
cx q37[22],q37[11];
rz(0.5853981633974479) q37[11];
sx q37[22];
cx q37[22],q37[11];
rz(-pi) q37[11];
u3(pi,-3*pi/4,3*pi/4) q37[22];
u3(pi,-3*pi/4,3*pi/4) q37[23];
sx q37[3];
cx q37[14],q37[3];
u2(0.9853981633974493,pi/2) q37[14];
rz(-0.9853981633974493) q37[3];
cx q37[14],q37[3];
sx q37[14];
rz(0.9853981633974496) q37[3];
cx q37[14],q37[3];
u2(pi/2,0) q37[14];
cx q37[15],q37[14];
rz(-0.8000000000000006) q37[14];
u2(0.8000000000000007,pi/2) q37[15];
cx q37[15],q37[14];
rz(0.8000000000000008) q37[14];
sx q37[15];
cx q37[15],q37[14];
u2(0,pi/2) q37[14];
u3(pi,-3*pi/4,3*pi/4) q37[15];
u2(-pi/2,-pi/2) q37[3];
cx q37[14],q37[3];
u2(0.5853981633974485,pi/2) q37[14];
rz(0.5853981633974487) q37[3];
cx q37[14],q37[3];
sx q37[14];
rz(0.5853981633974479) q37[3];
cx q37[14],q37[3];
u3(pi,-3*pi/4,3*pi/4) q37[14];
rz(-pi) q37[3];
sx q37[4];
u2(-pi,-pi/2) q37[5];
cx q37[5],q37[4];
rz(-0.9853981633974493) q37[4];
u2(0.9853981633974493,pi/2) q37[5];
cx q37[5],q37[4];
rz(0.9853981633974496) q37[4];
sx q37[5];
cx q37[5],q37[4];
u2(-pi/2,-pi/2) q37[4];
u2(pi/2,0) q37[5];
u2(-pi,-pi/2) q37[6];
cx q37[6],q37[5];
rz(-0.8000000000000006) q37[5];
u2(0.8000000000000007,pi/2) q37[6];
cx q37[6],q37[5];
rz(0.8000000000000008) q37[5];
sx q37[6];
cx q37[6],q37[5];
u2(0,pi/2) q37[5];
cx q37[5],q37[4];
rz(0.5853981633974487) q37[4];
u2(0.5853981633974485,pi/2) q37[5];
cx q37[5],q37[4];
rz(0.5853981633974479) q37[4];
sx q37[5];
cx q37[5],q37[4];
rz(-pi) q37[4];
u3(pi,-3*pi/4,3*pi/4) q37[5];
u3(pi,-3*pi/4,3*pi/4) q37[6];
sx q37[7];
cx q37[18],q37[7];
u2(0.9853981633974493,pi/2) q37[18];
rz(-0.9853981633974493) q37[7];
cx q37[18],q37[7];
sx q37[18];
rz(0.9853981633974496) q37[7];
cx q37[18],q37[7];
u2(pi/2,0) q37[18];
cx q37[19],q37[18];
rz(-0.8000000000000006) q37[18];
u2(0.8000000000000007,pi/2) q37[19];
cx q37[19],q37[18];
rz(0.8000000000000008) q37[18];
sx q37[19];
cx q37[19],q37[18];
u2(0,pi/2) q37[18];
u3(pi,-3*pi/4,3*pi/4) q37[19];
u2(-pi/2,-pi/2) q37[7];
cx q37[18],q37[7];
u2(0.5853981633974485,pi/2) q37[18];
rz(0.5853981633974487) q37[7];
cx q37[18],q37[7];
sx q37[18];
rz(0.5853981633974479) q37[7];
cx q37[18],q37[7];
u3(pi,-3*pi/4,3*pi/4) q37[18];
rz(-pi) q37[7];
sx q37[8];
u2(-pi,-pi/2) q37[9];
cx q37[9],q37[8];
rz(-0.9853981633974493) q37[8];
u2(0.9853981633974493,pi/2) q37[9];
cx q37[9],q37[8];
rz(0.9853981633974496) q37[8];
sx q37[9];
cx q37[9],q37[8];
u2(-pi/2,-pi/2) q37[8];
u2(pi/2,0) q37[9];
cx q37[10],q37[9];
u2(0.8000000000000007,pi/2) q37[10];
rz(-0.8000000000000006) q37[9];
cx q37[10],q37[9];
sx q37[10];
rz(0.8000000000000008) q37[9];
cx q37[10],q37[9];
u3(pi,-3*pi/4,3*pi/4) q37[10];
u2(0,pi/2) q37[9];
cx q37[9],q37[8];
rz(0.5853981633974487) q37[8];
u2(0.5853981633974485,pi/2) q37[9];
cx q37[9],q37[8];
rz(0.5853981633974479) q37[8];
sx q37[9];
cx q37[9],q37[8];
rz(-pi) q37[8];
u3(pi,-3*pi/4,3*pi/4) q37[9];
barrier q37[0],q37[1],q37[2],q37[3],q37[4],q37[5],q37[6],q37[7],q37[8],q37[9],q37[10],q37[11],q37[12],q37[13],q37[14],q37[15],q37[16],q37[17],q37[18],q37[19],q37[20],q37[21],q37[22],q37[23];
u2(-pi,-pi/2) q37[0];
sx q37[1];
sx q37[10];
u2(-pi,-pi/2) q37[11];
cx q37[11],q37[10];
rz(-0.9853981633974493) q37[10];
u2(0.9853981633974493,pi/2) q37[11];
cx q37[11],q37[10];
rz(0.9853981633974496) q37[10];
sx q37[11];
cx q37[11],q37[10];
u2(-pi/2,-pi/2) q37[10];
u2(pi/2,0) q37[11];
cx q37[0],q37[11];
u2(0.8000000000000007,pi/2) q37[0];
rz(-0.8000000000000006) q37[11];
cx q37[0],q37[11];
sx q37[0];
rz(0.8000000000000008) q37[11];
cx q37[0],q37[11];
u3(pi,-3*pi/4,3*pi/4) q37[0];
u2(0,pi/2) q37[11];
cx q37[11],q37[10];
rz(0.5853981633974487) q37[10];
u2(0.5853981633974485,pi/2) q37[11];
cx q37[11],q37[10];
rz(0.5853981633974479) q37[10];
sx q37[11];
cx q37[11],q37[10];
rz(-pi) q37[10];
u3(pi,-3*pi/4,3*pi/4) q37[11];
u2(-pi,-pi/2) q37[12];
cx q37[12],q37[1];
rz(-0.9853981633974493) q37[1];
u2(0.9853981633974493,pi/2) q37[12];
cx q37[12],q37[1];
rz(0.9853981633974496) q37[1];
sx q37[12];
cx q37[12],q37[1];
u2(-pi/2,-pi/2) q37[1];
u2(pi/2,0) q37[12];
u2(-pi,-pi/2) q37[13];
cx q37[13],q37[12];
rz(-0.8000000000000006) q37[12];
u2(0.8000000000000007,pi/2) q37[13];
cx q37[13],q37[12];
rz(0.8000000000000008) q37[12];
sx q37[13];
cx q37[13],q37[12];
u2(0,pi/2) q37[12];
cx q37[12],q37[1];
rz(0.5853981633974487) q37[1];
u2(0.5853981633974485,pi/2) q37[12];
cx q37[12],q37[1];
rz(0.5853981633974479) q37[1];
sx q37[12];
cx q37[12],q37[1];
rz(-pi) q37[1];
u3(pi,-3*pi/4,3*pi/4) q37[12];
u3(pi,-3*pi/4,3*pi/4) q37[13];
u2(-pi,-pi/2) q37[16];
u2(-pi,-pi/2) q37[17];
sx q37[2];
u2(-pi,-pi/2) q37[20];
u2(-pi,-pi/2) q37[21];
u2(-pi,-pi/2) q37[3];
cx q37[3],q37[2];
rz(-0.9853981633974493) q37[2];
u2(0.9853981633974493,pi/2) q37[3];
cx q37[3],q37[2];
rz(0.9853981633974496) q37[2];
sx q37[3];
cx q37[3],q37[2];
u2(-pi/2,-pi/2) q37[2];
u2(pi/2,0) q37[3];
u2(-pi,-pi/2) q37[4];
cx q37[4],q37[3];
rz(-0.8000000000000006) q37[3];
u2(0.8000000000000007,pi/2) q37[4];
cx q37[4],q37[3];
rz(0.8000000000000008) q37[3];
sx q37[4];
cx q37[4],q37[3];
u2(0,pi/2) q37[3];
cx q37[3],q37[2];
rz(0.5853981633974487) q37[2];
u2(0.5853981633974485,pi/2) q37[3];
cx q37[3],q37[2];
rz(0.5853981633974479) q37[2];
sx q37[3];
cx q37[3],q37[2];
rz(-pi) q37[2];
u3(pi,-3*pi/4,3*pi/4) q37[3];
u3(pi,-3*pi/4,3*pi/4) q37[4];
sx q37[5];
cx q37[16],q37[5];
u2(0.9853981633974493,pi/2) q37[16];
rz(-0.9853981633974493) q37[5];
cx q37[16],q37[5];
sx q37[16];
rz(0.9853981633974496) q37[5];
cx q37[16],q37[5];
u2(pi/2,0) q37[16];
cx q37[17],q37[16];
rz(-0.8000000000000006) q37[16];
u2(0.8000000000000007,pi/2) q37[17];
cx q37[17],q37[16];
rz(0.8000000000000008) q37[16];
sx q37[17];
cx q37[17],q37[16];
u2(0,pi/2) q37[16];
u3(pi,-3*pi/4,3*pi/4) q37[17];
u2(-pi/2,-pi/2) q37[5];
cx q37[16],q37[5];
u2(0.5853981633974485,pi/2) q37[16];
rz(0.5853981633974487) q37[5];
cx q37[16],q37[5];
sx q37[16];
rz(0.5853981633974479) q37[5];
cx q37[16],q37[5];
u3(pi,-3*pi/4,3*pi/4) q37[16];
rz(-pi) q37[5];
sx q37[6];
u2(-pi,-pi/2) q37[7];
cx q37[7],q37[6];
rz(-0.9853981633974493) q37[6];
u2(0.9853981633974493,pi/2) q37[7];
cx q37[7],q37[6];
rz(0.9853981633974496) q37[6];
sx q37[7];
cx q37[7],q37[6];
u2(-pi/2,-pi/2) q37[6];
u2(pi/2,0) q37[7];
u2(-pi,-pi/2) q37[8];
cx q37[8],q37[7];
rz(-0.8000000000000006) q37[7];
u2(0.8000000000000007,pi/2) q37[8];
cx q37[8],q37[7];
rz(0.8000000000000008) q37[7];
sx q37[8];
cx q37[8],q37[7];
u2(0,pi/2) q37[7];
cx q37[7],q37[6];
rz(0.5853981633974487) q37[6];
u2(0.5853981633974485,pi/2) q37[7];
cx q37[7],q37[6];
rz(0.5853981633974479) q37[6];
sx q37[7];
cx q37[7],q37[6];
rz(-pi) q37[6];
u3(pi,-3*pi/4,3*pi/4) q37[7];
u3(pi,-3*pi/4,3*pi/4) q37[8];
sx q37[9];
cx q37[20],q37[9];
u2(0.9853981633974493,pi/2) q37[20];
rz(-0.9853981633974493) q37[9];
cx q37[20],q37[9];
sx q37[20];
rz(0.9853981633974496) q37[9];
cx q37[20],q37[9];
u2(pi/2,0) q37[20];
cx q37[21],q37[20];
rz(-0.8000000000000006) q37[20];
u2(0.8000000000000007,pi/2) q37[21];
cx q37[21],q37[20];
rz(0.8000000000000008) q37[20];
sx q37[21];
cx q37[21],q37[20];
u2(0,pi/2) q37[20];
u3(pi,-3*pi/4,3*pi/4) q37[21];
u2(-pi/2,-pi/2) q37[9];
cx q37[20],q37[9];
u2(0.5853981633974485,pi/2) q37[20];
rz(0.5853981633974487) q37[9];
cx q37[20],q37[9];
sx q37[20];
rz(0.5853981633974479) q37[9];
cx q37[20],q37[9];
u3(pi,-3*pi/4,3*pi/4) q37[20];
rz(-pi) q37[9];
barrier q37[0],q37[1],q37[2],q37[3],q37[4],q37[5],q37[6],q37[7],q37[8],q37[9],q37[10],q37[11],q37[12],q37[13],q37[14],q37[15],q37[16],q37[17],q37[18],q37[19],q37[20],q37[21],q37[22],q37[23];