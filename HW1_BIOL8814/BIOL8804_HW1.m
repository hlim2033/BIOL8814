% C Cultures and each has N cells. Mutation prob = mu. 

%% PROBLEM 1 : ONE GENERATION OF LD EXPERIMENT
% C = 500, N = 1000, mu = 1*10^-3
clear
clc 

C = 500;
N = 1000;
mu = 1*10^-3;

for i = 1:C
    i;
    ct = 0;
    for j = 1:N
        for k = 1:2
            ct = ct + 1; 
            ALREADYONE = 0;
            rn = rand;
            if ALREADYONE == 1
                continue
            end

            if rn < mu
                state = 1;
                ALREADYONE = 1;
            else
                state = 0;
            end
            STATE(ct) = state;
        end
    end
    SUM = sum(STATE);
    SIMUL(i) = SUM;
end
figure
hist(SIMUL)

No_Zeros = sum(SIMUL == 0) 
lambda = whatislambda(SIMUL) %labmda = mu * 2 * N (after dividing) 
mean(SIMUL)
var(SIMUL)
%%
[y, x] = hist(SIMUL)
plot(y)

%% PROBLEM 2 
% C = 1000, N = 400, mu = 1*10^-7
clear
clc 
C = 1000;
N = 400;
mu = 1*10^-7;
i;
C=1000;
for i = 1:C
    i
    Original_STATE = zeros(1,400);
    age_ct = 0;
    for g = 1:15
        ct = 0;
        for j = 1:length(Original_STATE)
            cell = Original_STATE(j);
            if cell == 0 
                for k = 1:2
                    ct = ct + 1;
                    ALREADYONE = 0;
                    rn = rand;
                    if ALREADYONE == 1
                        continue
                    end

                    if rn < mu
                        state = 1;
                        ALREADYONE = 1;
                    else
                        state = 0;
                    end
                    STATE(ct) = state;
                end
            elseif cell == 1
                for q = 1:2
                    ct = ct + 1;
                    state = 1;
                    STATE(ct) = state;
                end
            end
        end
        Original_STATE = STATE;
        ct;
        if sum(STATE) > 0 && age_ct ==0
            age = g; 
            age_ct = age_ct+1;
        end
    end
    SUM = sum(Original_STATE);
    SIMUL(i) = SUM;
    AGE(i) = age;
    clear STATE
    clear Original_STATE
end
SIMUL
AGE
%% Problem 2, CONT'D
clear
clc
load prob2_variable.mat
DATA = SIMUL;
DATA(287) = [];
n_D = length(DATA)

hist(DATA,1:2054)

DATA_sorted = sort(DATA)
no_jp_DATA = DATA_sorted(1:round(n_D*0.95))
var(no_jp_DATA)

%% PROBLEM3 AGE
clear
clc
load prob2_variable.mat
DATA = SIMUL;
DATA(287) = [];
age_DATA = AGE;
age_DATA(287) = [];
figure
hold on 
plot(DATA)
plot(age_DATA)
ylim([0,100])
hold off 

figure
scatter(age_DATA, DATA)

%% Problem 4
load prob2_variable.mat
p4_DATA = SIMUL(1:100);
p4_AGE = AGE(1:100);

figure
hist(p4_DATA)

figure
hist(DATA, 150)
lambda = whatislambda(DATA) %Poission? NO!

No_Zeros = sum(DATA == 0) 
%No_nonZeros = sum(DATA ~= 0) 

mean_mutants = mean(DATA)