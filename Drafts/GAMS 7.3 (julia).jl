,
using JuMP
using Ipopt
using GLPKMathProgInterface
using Cbc
using DataArrays, DataFrames, CSV
using Gadfly
using Pajarito#using CPLEX

###  Data definition
#Generators data

GenData = DataFrame(
    Bus = [18,21,1,2,15,16,23,23,7,13,15,22],
    Pmax = [400,400,152,152,155,155,310,350,350,591,60,300],
    Pmin = [100,100,30.40,30.40,54.25,54.25,108.50,140,75,206.85,12,0],
    c = [5.47,5.47,13.32,13.32,16,10.52,10.52,10.89,20.7,20.93,26.11,0],
    a = [0.5,0.6,0.7,0.4,0.5,0.7,0.5,0.3,0.4,0.5,0.7,0.8],
    RU = [47,47,14,14,21,21,21,28,49,21,7,35],
    RD = [47,47,14,14,21,21,21,28,49,21,7,35])

#Number of generators
NGen = size(GenData,1)
# Generator limits
Pmin = zeros(NGen);
Pmax = zeros(NGen);
Qmin = zeros(NGen);
Qmax = zeros(NGen);
# Generation costs
c = zeros(NGen);
a = zeros(NGen);
Bus = zeros(NGen);
#RumpUp
RU = zeros(NGen);
#RumpDown
RD = zeros(NGen);



for i=1:NGen
    Bus[i] = GenData[i,:Bus]
    Pmin[i] = GenData[i,:Pmin]
    Pmax[i] = GenData[i,:Pmax]
    c[i] = GenData[i,:c]
    a[i] = GenData[i,:c]
    RU[i] = GenData[i,:RU]
    RD[i] = GenData[i,:RD]

end

# Lines data
BranchData = DataFrame(
    FromBus = [1,1,1,2,2,3,3,4,5,6,7,8,8,9,9,10,10,11,11,12,12,13,14,15,15,15,16,16,17,17,18,19,20,21],
    ToBus = [2,3,5,4,6,9,24,9,10,10,8,9,10,11,12,11,12,13,14,13,23,23,16,16,21,24,17,19,18,22,21,20,23,22],
    r = [0.0026,0.0546,0.0218,0.0328,0.0497,0.0308,0.0023,0.0268,0.0228,0.0139,0.0159,0.0427,0.0427,0.0023,0.0023,0.0023,0.0023,0.0061,0.0054,0.0061,0.0124,0.0111,0.0050,0.0022,0.0032,0.0067,0.0033,0.0030,0.0018,0.0135,0.0017,0.0026,0.0014,0.0087],
    x = [0.0139,0.2112,0.0845,0.1267,0.1920,0.1190,0.0839,0.1037,0.0883,0.0605,0.0614,0.1651,0.1651,0.0839,0.0839,0.0839,0.0839,0.0476,0.0418,0.0476,0.0966,0.0865,0.0389,0.0173,0.0245,0.0519,0.0259,0.0231,0.0144,0.1053,0.0130,0.0198,0.0108,0.0678],
    b = [0.4611,0.0572,0.0229,0.0343,0.0520,0.0322,0.0000,0.0281,0.0239,2.4590,0.0166,0.0447,0.0447,0.0000,0.0000,0.0000,0.0000,0.0999,0.0879,0.0999,0.2030,0.1818,0.0818,0.0364,0.2060,0.1091,0.0545,0.0485,0.0303,0.2212,0.1090,0.1666,0.0910,0.1424],
    Smax = [175,175,175,175,175,175,400,175,175,175,175,175,175,400,400,400,400,500,500,500,500,500,500,500,1000,500,500,500,500,500,1000,1000,1000,500],)


NBuses = max(maximum(BranchData[:,:FromBus]), maximum(BranchData[:,:ToBus]));
NLines = size(BranchData,1);
BranchData[:Z]=BranchData[:r]+im*BranchData[:x];
BranchData[:Y]=1./BranchData[:Z];
Ybus = zeros(Complex,NBuses,NBuses);
Nodes = zeros(Int8,NLines,2); #node's matrix
SLmax = zeros(NLines); #line limits
x = zeros(NLines);
r = zeros(NLines);


for i=1:NLines

    Ybus[BranchData[i,:FromBus],BranchData[i,:ToBus]]=-BranchData[i,:Y]
    Ybus[BranchData[i,:ToBus],BranchData[i,:FromBus]]=-BranchData[i,:Y]
    Nodes[i,1] = BranchData[i,:FromBus]
    Nodes[i,2] = BranchData[i,:ToBus]
    SLmax[i] = BranchData[i,:Smax]
    x[i] = BranchData[i,:x]
    r[i] = BranchData[i,:x]
end

for i=1:NBuses
    for j=1:NLines
        if (i==BranchData[j,:FromBus])||(i==BranchData[j,:ToBus])
            Ybus[i,i]+=BranchData[j,:Y]
        end
    end
end

Gbus = real(Ybus);
Bbus = imag(Ybus);

NodesData = DataFrame(
    Pd = [108,97,180,74,71,136,125,171,175,195,0,0,265,194,317,100,0,333,181,128,0,0,0,0],
    Qd = [22,20,37,15,14,28,25,35,36,40,0,0,54,39,64,20,0,68,37,26,0,0,0,0],
    Wind = [0,0,0,0,0,0,0,200,0,0,0,0,0,0,0,0,0,0,150,0,100,0,0,0],
    SOC = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,200,0,100,0,0,0],
)
ScaleData = DataFrame(
    w = [0.0786666666666667,0.0866666666666667,0.117333333333333,0.258666666666667,0.361333333333333,0.566666666666667,0.650666666666667,0.566666666666667,0.484,0.548,0.757333333333333,0.710666666666667,0.870666666666667,0.932,0.966666666666667,1,0.869333333333333,0.665333333333333,0.656 ,0.561333333333333,0.565333333333333,0.556,0.724,0.84],
    d = [0.684511335492475,0.644122690036197,0.61306915602972,0.599733282530006,0.588874071251667,0.5980186702229,0.626786054486569,0.651743189178891,0.706039245570585,0.787007048961707,0.839016955610593,0.852733854067441,0.870642027052772,0.834254143646409,0.816536483139646,0.819394170318156,0.874071251666984,1,0.983615926843208,0.936368832158506,0.887597637645266,0.809297008954087,0.74585635359116,0.733473042484283],
)

T = size(ScaleData,1); #Time scale

Pload = zeros(NBuses); # Demand
Ql = zeros(NBuses);
Wind = zeros(NBuses); #Wind turbaine capacity
SOC_d = zeros(NBuses); #SOC capacity

w = zeros(T); #Wind profile coefficients
d = zeros(T); #Demand profile coefficients


for i=1:NBuses

    Pload[i] = NodesData[i,:Pd]
    Ql[i] = NodesData[i,:Qd]
    Wind[i] = NodesData[i,:Wind]
    SOC_d[i] = NodesData[i,:SOC]

end

for i=1:T
    w[i] = ScaleData[i,:w]
    d[i] = ScaleData[i,:d]
end

# Voltage limits
Vmin = 0.9;
Vmax = 1.15;

#Angle limits

θmax = pi/2;
θmin = -pi/2;

SlackBus = 13;

Sbase = 100;

SOC_max = 20;
SOC_min = 0.2*SOC_max;
SOC_0 = 0.2*SOC_max; #initial charge
Pd_min = 0*SOC_max; #discharge limits
Pd_max = 0.2*SOC_max;
Pc_min = 0*SOC_max; #charge limits
Pc_max = 0.2*SOC_max;
#charge/discharge efficiency
ηс = 0.95;
ηd = 0.9;
#Wind power limits
Pw_min = 0*Wind;
Pw_max = Wind;
Pwc_min = 0*Wind;
Pwc_max = Wind;
#Load shedding limits
Ls_min = 0*Pload;
Ls_max = Pload;

VOLL = 10000;
VOLW = 50;

N_max = 15; #total number of ESS in the system
NESS_max = 5; #limit of ESS in the bus

# Set size
NI = NBuses;
NL = NLines;
NG = NGen;

m = Model(solver = GLPKSolverMIP());
#m = Model(solver = IpoptSolver());
#m = Model(solver =  CbcSolver());
#m = Model(solver = PajaritoSolver(mip_solver = GLPKSolverMIP(), cont_solver = IpoptSolver()))

@variable(m, p[g=1:NG,t=1:T]);
@variable(m, pl[l=1:NL,t=1:T]);
@variable(m, θ[i=1:NI,t=1:T]);
@variable(m, Ls[i=1:NI,t=1:T]);
@variable(m, SOC[i=1:NI,t=1:T]);
@variable(m, pd[i=1:NI,t=1:T]);
@variable(m, pc[i=1:NI,t=1:T]);
@variable(m, pw[i=1:NI,t=1:T]);
@variable(m, pwc[i=1:NI,t=1:T]);
@variable(m, NESS[i=1:NI], Int);
#@variable(m, NESS[i=1:NI]);

#Nodal balance equation
@constraint(m, BalanceEq[i=1:NI, t=1:T], sum(p[g,t] for g=1:NG if Bus[g]==i) + Ls[i,t] + pw[i,t] - Pload[i]*d[t]/Sbase - pc[i,t] + pd[i,t] == sum(pl[l,t] for l=1:NL if Nodes[l,1]==i) - sum(pl[l,t] for l=1:NL if Nodes[l,2]==i));

#Branch power flow defenition
@constraint(m, PowerFlow[l=1:NL,t=1:T], pl[l,t] == (θ[Nodes[l,1],t]-θ[Nodes[l,2],t])/x[l]);
@constraint(m, FlowLimits[l=1:NL,t=1:T], -SLmax[l]/Sbase <= pl[l,t] <= SLmax[l]/Sbase);

#Technical generation limits
@constraint(m, GenLimitsP[g=1:NG,t=1:T], Pmin[g]/Sbase <= p[g,t] <= Pmax[g]/Sbase);
@constraint(m, RumpUp[g=1:NG,t=2:T], p[g,t] - p[g,t-1]  <= RU[g]/Sbase);
@constraint(m, RumpDown[g=1:NG,t=2:T], p[g,t-1] - p[g,t]  <= RD[g]/Sbase);

#Angle limits
@constraint(m, limitAngle1[i=SlackBus,t=1:T], θ[i,t] == 0);
@constraint(m, limitAngle2[i=1:NI,t=1:T], θmin <= θ[i,t] <= θmax);

#Charge balance
@constraint(m,ChargeBal1[i=1:NI,t=1], SOC[i,t] == NESS[i]*SOC_0/Sbase  + pc[i,t]*ηс - pd[i,t]/ηd);
@constraint(m,ChargeBal2[i=1:NI,t=2:T], SOC[i,t] == SOC[i,t-1] + pc[i,t]*ηс - pd[i,t]/ηd);
@constraint(m,ChargeBal3[i=1:NI,t=T], SOC[i,t] == NESS[i]*SOC_0/Sbase);
@constraint(m,ChargeBal4[i=1:NI,t=1:T], SOC[i,t] <= NESS[i]*SOC_max/Sbase);
@constraint(m,ChargeBal5[i=1:NI,t=1:T], pd[i,t] <= NESS[i]*Pd_max/Sbase);
@constraint(m,ChargeBal6[i=1:NI,t=1:T], pc[i,t] <= NESS[i]*Pc_max/Sbase);
@constraint(m,ChargeBal7[i=1:NI,t=1:T], SOC[i,t] >= NESS[i]*SOC_min/Sbase);
@constraint(m,ChargeBal8[i=1:NI,t=1:T], pd[i,t] >= NESS[i]*Pd_min/Sbase);
@constraint(m,ChargeBal9[i=1:NI,t=1:T], pc[i,t] >= NESS[i]*Pc_min/Sbase);
@constraint(m,ChargeBal10, sum(NESS[i] for i=1:NI) <= N_max);
@constraint(m,ChargeBal11[i=1:NI], NESS[i]  <= NESS_max);

#Wind constraine
@constraint(m,Wind1[i=1:NI,t=1:T], pw[i,t] + pwc[i,t] == w[t]*Wind[i]/Sbase);
@constraint(m,Wind2[i=1:NI,t=1:T], Pw_min[i]/Sbase <= pw[i,t] <= w[t]*Pw_max[i]/Sbase);
@constraint(m,Wind3[i=1:NI,t=1:T], Pwc_min[i]/Sbase <= pwc[i,t] <= w[t]*Pwc_max[i]/Sbase);

#Load shedding
@constraint(m,Load[i=1:NI,t=1:T], Ls_min[i]/Sbase <= Ls[i,t] <= d[t]*Ls_max[i]/Sbase);

#Objective function
@objective(m, Min, sum(c[g]*p[g,t]*Sbase for g=1:NG,t=1:T) + sum(VOLL*Ls[i,t]*Sbase for i=1:NI,t=1:T) + sum(VOLW*pwc[i,t]*Sbase for i=1:NI,t=1:T));

solve(m)

OF = getobjectivevalue(m)

NEES_f = getvalue(NESS)

Pg=getvalue(p);
GenPower = DataFrame();
GenPower[:g1] = Pg[1,:];
GenPower[:g2] = Pg[2,:];
GenPower[:g3] = Pg[3,:];
GenPower[:g4] = Pg[4,:];
GenPower[:g5] = Pg[5,:];
GenPower[:g6] = Pg[6,:];
GenPower[:g7] = Pg[7,:];
GenPower[:g8] = Pg[8,:];
GenPower[:g9] = Pg[9,:];
GenPower[:g10] = Pg[10,:];
GenPower[:g11] = Pg[11,:];
GenPower[:g12] = Pg[12,:]

plot(GenPower, x = Row.index, y=Col.value(:g1,:g2,:g3,:g4,:g5,:g6,:g7,:g8,:g9,:g10,:g11,:g12),
    color=Col.index(:g1,:g2,:g3,:g4,:g5,:g6,:g7,:g8,:g9,:g10,:g11,:g12),
    Geom.line,

Guide.XLabel("Time"),
Guide.YLabel("Power"),
#Guide.Title("Title")
)

PowerBalance = DataFrame();
PowerBalance[:Load] = sum(Pload)*d[:];
PowerBalance[:Gen] = sum(getvalue(p[g,:]) for g=1:NG)*Sbase;
PowerBalance[:SOC] = sum(getvalue(SOC[i,:]) for i=1:NI)*Sbase;
PowerBalance[:Wind] = sum(getvalue(pw[i,:]) for i=1:NI)*Sbase;

plot(PowerBalance, x = Row.index, y=Col.value(:Load,:Gen,:SOC,:Wind),
    color=Col.index(:Load,:Gen,:SOC,:Wind),
    Geom.line,

Guide.XLabel("Time"),
Guide.YLabel("Power"),
#Guide.Title("Title")
)
