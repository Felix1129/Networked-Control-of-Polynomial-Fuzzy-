 
clear all
clc
tic

count_1=0;%�p�ƥΡA�����ĴX�����
count_2=0;%�p�ƥΡA�����X��=0���
count_3=0;%�p�ƥΡA�����X��=0���
%�x�}�x�s��ơA�ݭn���ŧi�x�}�j�p

scalar1=12;
scalar2=12;
scalar3=12;
scalar4=12;
scalar5=12;
scalar_all = scalar1*scalar2*scalar3*scalar4*scalar5;
global fin_K1_data  fin_K2_data sol_space

%�q(0,0,0,0),(0,0,0,0.25),(0,0,0,0.5),(0,0,0,0.75)��(1,1,1,1)
%���F�N�Ҧ�gain���x�s�A�ŧifin_K_data��(scalar_all*2,6)��0�x�}
fin_K1_data=zeros(scalar_all*2,6);
fin_K2_data=zeros(scalar_all*2,6);
%�P�ˡA�Nm,n,lambda�s�isol_space�A�᭱�Q�α��󦡧P�O�A
%��S�x�Ȥj��0�A�x�s5�ѼƭȡA�H���e�X�ѪŶ�
sol_space=zeros(scalar_all,5);


%�ѩ���q���ȴ��դҦ��e�P�e���ܰѼ�lambda_s����0~1����
num = 5;                            % �ŧi�@�H���`�ƭ�

%%
for m=0.1:0.1:1
    for n=0.1:0.1:1
        for lambda=0.02:0.02:0.2
            for gamma=0.1:0.1:1
                for kappa=1:1:10
                    try
                        %% �L�X�����ܶq
                        for block=8
                            m
                            n
                            lambda
                            gamma
                            kappa
                            N=2
                        end
                        %% �ŧi�ܼƲŸ�
                        for block=1
                            % �ŧi�W���ܼ�
                            syms x1 x2 x3 x4 x5 x6
                            syms v11 v12 v13 v14 v15 v16 v17 v18 v19 v110 v111 v112
                            syms v21 v22 v23 v24 v25 v26 v27 v28 v29 v210 v211 v212 v213 v214 v215 v216 v217 v218 v219 v220 v221
                            syms v222 v223 v224  v225 v226 v227 v228 v229 v230 v231 v232 v233 v234 v235 v236 v237 v238
                            vars = [x1;x2;x3;x4;x5;x6];
                            s1=[v11;v12;v13;v14;v15;v16];
                            s2=[v17;v18;v19;v110;v111;v112];
                            %�]��phi�O15*15
                            s3=[v21; v22; v23; v24; v25; v26; v27; v28; v29; v210 ;v211 ;v212;v213;v214;v215;v216;v217;v218;v219;v220;v221;
                                v222; v223; v224; v225; v226; v227; v228; v229; v230; v231; v232; v233; v234; v235; v236; v237; v238];
                            x_h=[x1;x2;x3;x4;x5;x6];
                            x=x_h;
                            % �ŧi�t�ίx�}
                            A1 = [0  0  1  0  0  0;         % �ŧi�t�ίx�}
                                0  0  0  1  0  0;
                                0  0  0  0  0  0;
                                0  0  0  0  0  0;
                                1  0  0  0  0  0;
                                0  1  0  0  0  0];
                            A2 = A1;
                            
                            B1 = [0     0;                  % �ŧi��J�x�}
                                0     0;
                                4.8     0;
                                0   4.8;
                                0     0;
                                0     0];
                            B2 = [0     0;
                                0     0;
                                14.8     0;
                                0  14.8;
                                0     0;
                                0     0];
                            
                            C1 = [1  1  0  0  0  0];        % �ŧi��X�x�}
                            C2 = C1;
                            
                            D1=[0.5 ; 0.5 ; 0 ; 0; 0; 0];
                            D2=D1;
                            T=eye(6);%eye�O���x�},(n)���Ʀr�On*n���N��
                            deltaA1=[  0    0    1    0  0  0;
                                0    0    0    1  0  0;
                                0    0    0    0  0  0;
                                0    0    0    0  0  0;
                                1    0    0    0  0  0;
                                0    1    0    0  0  0];
                            deltaA2=[  0    0    1    0  0  0;
                                0    0    0    1  0  0;
                                0    0    0    0  0  0;
                                0    0    0    0  0  0;
                                1    0    0    0  0  0;
                                0    1    0    0  0  0];
                            deltaB1=[0     0;
                                0     0;
                                0     0;
                                0     0;
                                0     0;
                                0     0];
                            
                            deltaB2=deltaB1;
                            J=eye(6);
                            K=sin(90);
                            
                            Ra1=inv(K)*inv(J)*deltaA1;
                            Ra2=inv(K)*inv(J)*deltaA2;
                            Rb1=inv(K)*inv(J)*deltaB1;
                            Rb2=inv(K)*inv(J)*deltaB2;
                            
                            Raij=m*Ra1+n*Ra2;
                            Rbij=m*Rb1+n*Rb2;
                        end
                        %% Output Feedback    %�]����X�^�¡A�ݭn�ϥΥ���ɶi���ܼ��ܴ�
                        for block=2
                            L1 = [C1'*(inv(C1*C1')) ortc(C1')];     % �ŧi��X�^�¤��D�_���ܴ��x�}
                            L2 = [C2'*(inv(C2*C2')) ortc(C2')];
                            TA1 = L1\A1*L1;                           % �ŧi�ܴ��ᤧ�t�ίx�}
                            TA2 = L2\A2*L2;
                            
                            TB1 = L1\B1;                              % �ŧi�ܴ��ᤧ��J�x�}
                            TB2 = L2\B2;
                            
                            L1=eye(6);
                            L2=eye(6);
                            sizeA =[6 6];                                                                  % �ŧi�t�ίx�}������
                            sizeB = [2 6];                                                                 % �ŧi��J�x�}������
                            eps1=0.01;
                            eps2=0.01;
                            eps3=0.01;
                            J_tear=m*inv(L1)*J+n*inv(L1)*J;
                        end
                        %% �Ыئh�����x�}
                        for block=3 %��l�ƿW���ܼ�vars1, s1, s2, s3
                            prog=sosprogram([vars; s1; s2; s3]);
                            % �Ы�P1�x�}�A�䤤���|�ӿ�J�A�e�T�ӥ��n��J�H�γ̫�@�ӥi���J
                            %�Ĥ@�ӿ�J��SOS�{���A�ĤG�Ӭ��涵���V�q�A�ĤT�Ӭ�P�x�}�����סA�̫�@�Ӭ����wP�O�_����ٯx�}
                            
                            
                            %N=1�ݭn���x�}
                            [prog, P1] = sospolymatrixvar(prog,monomials(vars,0),sizeA,'symmetric');            % �Ы�P1�x�}
                            [prog, X1] = sospolymatrixvar(prog,monomials(vars,0),sizeA,'symmetric');            % �Ы�X1�x�}
                            [prog, Q1] = sospolymatrixvar(prog,monomials(vars,0),sizeA,'symmetric');           % �Ы�Q1�x�}
                            [prog,H11] = sospolymatrixvar(prog,monomials(vars,0),sizeA);                                 % �Ы�H11�x�}
                            [prog,H21] = sospolymatrixvar(prog,monomials(vars,0),sizeA);                                 % �Ы�H21�x�}
                            [prog,H31] = sospolymatrixvar(prog,monomials(vars,0),sizeA);                                 % �Ы�H31�x�}
                            [prog,M11] = sospolymatrixvar(prog,monomials(vars,0),sizeB);                                % �Ы�M11�x�}
                            [prog,M21] = sospolymatrixvar(prog,monomials(vars,0),sizeB);                                % �Ы�M21�x�}
                            %N=2�ݭn���x�}
                            %N=2�ɡA�W�U���x�}���ݭn
                            [prog, P2] = sospolymatrixvar(prog,monomials(vars,0),sizeA,'symmetric');            % �Ы�P2�x�}
                            [prog, Q2] = sospolymatrixvar(prog,monomials(vars,0),sizeA,'symmetric');           % �Ы�Q2�x�}
                            [prog, X2] = sospolymatrixvar(prog,monomials(vars,0),sizeA,'symmetric');            % �Ы�X2�x�}
                            [prog,H12] = sospolymatrixvar(prog,monomials(vars,0),sizeA);                                 % �Ы�H12�x�}
                            [prog,H22] = sospolymatrixvar(prog,monomials(vars,0),sizeA);                                 % �Ы�H22�x�}
                            [prog,H32] = sospolymatrixvar(prog,monomials(vars,0),sizeA);                                 % �Ы�H32�x�}
                            [prog,M12] = sospolymatrixvar(prog,monomials(vars,0),sizeB);                                 % �Ы�M12�x�}
                            [prog,M22] = sospolymatrixvar(prog,monomials(vars,0),sizeB);                                 % �Ы�M22�x�}
                        end
                        %% í�w����
                        for block=4
                            % �bSOS�{���[�J�x�}�����������A�������P1-eps1�x�}������
                            %eye(sizeA)�O�]��eps�O�`�ơA�x�}�����`��
                            %���w�zí�w����Ĥ@��
                            SOS1 = s1.'*(P1-eps1.*eye(sizeA))*s1;
                            prog = sosmatrixineq(prog,SOS1);
                            % �bSOS�{���[�J�x�}�����������A�������Q1-eps2�x�}������
                            %���w�zí�w����ĤG��
                            SOS2 = s2.'*(Q1-eps2.*eye(sizeA))*s2;
                            prog = sosmatrixineq(prog,SOS2);
                            %N=1
                            Sum11 = 0;                                                              % �ŧiP1�x�}��ɶ��L��(i=1)���M
                            Sum21 = 0;                                                              % �ŧiP1�x�}��ɶ��L��(i=2)���M
                            
                            for i = 1:length(vars)                                                  % �̾�x_tilde���ƶq(i=1)�i��X��for�j��
                                Sum11 = Sum11 + diff(P1,vars(i))*(TA1(i,:)*vars(1:sizeA));          % �p��P1�x�}��ɶ��L��(i=1)�ìۥ[
                            end
                            
                            for i = 1:length(vars)                                                  % �̾�x_tilde���ƶq(i=2)�i��X��for�j��
                                Sum21 = Sum21 + diff(P1,vars(i))*(TA2(i,:)*vars(1:sizeA));          % �p��P1�x�}��ɶ��L��(i=2)�ìۥ[
                            end
                            %11111
                            %
                            for subblock=1
                                %% í�w����3-52,i=j=g=s=1,(4*3)x(4*3):
                                %�]��i=j=g=s=1
                                % Aij_tear=inv(L1)*T*(mAi+nAj)*L1;
                                % Bij_tear=inv(L1)*T*(mAi+nAj);
                                % Gij_tear=mCi+nCj;
                                % Dij_tear=inv(L1)*T*(mDi+nDj);
                                % Mgse_tear=mMge+nMse;
                                Aij_tear=inv(L1)*T*(m*A1+n*A1)*L1;
                                Bij_tear=inv(L1)*T*(m*B1+n*B1);
                                Cij_tear=(m*C1+n*C1)*L1;
                                Dij_tear=inv(L1)*T*(m*D1+n*D1);
                                Mgse_tear=m*M11+n*M11;
                                Sum=Sum11;
                                
                                Oa111=diff(P1,x2).*(A1(3,:)*x_h)+H11+H11.'-Aij_tear*X1-X1.'*Aij_tear+Sum;
                                Oa121=-H11+H21.'-Bij_tear*Mgse_tear-X1.'*Aij_tear.';
                                Oa131=P1+H31.'+X1-X1.'*Aij_tear;
                                Oa141=-Dij_tear;
                                Oa151=lambda*H11;
                                Oa161=(Cij_tear*X1).';
                                Oa171=-J_tear*X1;
                                Oa181=kappa*X1.'*Raij.';
                                
                                Oa211=Oa121.';
                                Oa221=-H21-H21.'-Bij_tear*Mgse_tear-Mgse_tear.'*Bij_tear.';
                                Oa231=-H31.'+X1-Mgse_tear.'*Bij_tear.';
                                Oa241=-Dij_tear;
                                Oa251=lambda*H21;
                                Oa261=zeros(6,1);
                                Oa271=-J_tear*X1;
                                Oa281=kappa* Mgse_tear.'*Rbij.';
                                
                                Oa311=Oa131.';
                                Oa321=Oa231.';
                                Oa331=lambda*Q1+X1.'+X1;
                                Oa341=-Dij_tear;
                                Oa351=lambda*H31;
                                Oa361=zeros(6,1);
                                Oa371=-J_tear*X1;
                                Oa381=zeros(6,6);
                                
                                Oa411=Oa141.';
                                Oa421=Oa241.';
                                Oa431=Oa341.';
                                Oa441=-gamma*gamma*eye(1);
                                Oa451=zeros(1,6);
                                Oa461=zeros(1,1);
                                Oa471=zeros(1,6);
                                Oa481=zeros(1,6);
                                
                                Oa511=Oa151.';
                                Oa521=Oa251.';
                                Oa531=Oa351.';
                                Oa541=Oa451.';
                                Oa551=-lambda*Q1;
                                Oa561=zeros(6,1);
                                Oa571=zeros(6,6);
                                Oa581=zeros(6,6);
                                
                                
                                
                                Oa611=Oa161.';
                                Oa621=Oa261.';
                                Oa631=Oa361.';
                                Oa641=Oa461.';
                                Oa651=Oa561.';
                                Oa661=-eye(1);
                                Oa671=zeros(1,6);
                                Oa681=zeros(1,6);
                                
                                Oa711=Oa171.';
                                Oa721=Oa271.';
                                Oa731=Oa371.';
                                Oa741=Oa471.';
                                Oa751=Oa571.';
                                Oa761=Oa671.';
                                Oa771=-kappa*eye(6);
                                Oa781=zeros(6,6);
                                
                                Oa811=Oa181.';
                                Oa821=Oa281.';
                                Oa831=Oa381.';
                                Oa841=Oa481.';
                                Oa851=Oa581.';
                                Oa861=Oa681.';
                                Oa871=Oa781.';
                                Oa881=-kappa*eye(6);
                                
                                
                                
                                O11111=[Oa111  Oa121  Oa131  Oa141  Oa151  Oa161  Oa171  Oa181;
                                    Oa211  Oa221  Oa231  Oa241  Oa251  Oa261  Oa271  Oa281;
                                    Oa311  Oa321  Oa331  Oa341  Oa351  Oa361  Oa371  Oa381;
                                    Oa411  Oa421  Oa431  Oa441  Oa451  Oa461  Oa471  Oa481;
                                    Oa511  Oa521  Oa531  Oa541  Oa551  Oa561  Oa571  Oa581;
                                    Oa611  Oa621  Oa631  Oa641  Oa651  Oa661  Oa671  Oa681;
                                    Oa711  Oa721  Oa731  Oa741  Oa751  Oa761  Oa771  Oa781;
                                    Oa811  Oa821  Oa831  Oa841  Oa851  Oa861  Oa871  Oa881];
                                SOS3=-s3.'*([O11111]+eps3.*eye(38))*s3;
                                %                             SOS3=-s3.'*([O11111]+eps3.*eye(38))*s3;
                                % �ŧi�������x�}
                                
                                prog = sosineq(prog,SOS3,'sparsemultipartite',{ vars.' s3.'});
                            end
                            %}
                            
                            
                            %22221
                            for subblock=2
                                %% í�w����3-52,i=j=g=s=2,(4*3)x(4*3):
                                %�]��i=j=g=s=2
                                % Aij_tear=inv(L1)*T*(mAi+nAj)*L1;
                                % Bij_tear=inv(L1)*T*(mAi+nAj);
                                % Gij_tear=mCi+nCj;
                                % Dij_tear=inv(L1)*T*(mDi+nDj);
                                % Mgse_tear=mMge+nMse;
                                Aij_tear=inv(L1)*T*(m*A2+n*A2)*L1;
                                Bij_tear=inv(L1)*T*(m*B2+n*B2);
                                Cij_tear=(m*C2+n*C2)*L1;
                                Dij_tear=inv(L1)*T*(m*D2+n*D2);
                                Mgse_tear=m*M21+n*M21;
                                Sum=Sum21;
                                
                                Oa111=diff(P1,x2).*(A1(3,:)*x_h)+H11+H11.'-Aij_tear*X1-X1.'*Aij_tear+Sum;
                                Oa121=-H11+H21.'-Bij_tear*Mgse_tear-X1.'*Aij_tear.';
                                Oa131=P1+H31.'+X1-X1.'*Aij_tear;
                                Oa141=-Dij_tear;
                                Oa151=lambda*H11;
                                Oa161=(Cij_tear*X1).';
                                Oa171=-J_tear*X1;
                                Oa181=kappa*X1.'*Raij.';
                                
                                Oa211=Oa121.';
                                Oa221=-H21-H21.'-Bij_tear*Mgse_tear-Mgse_tear.'*Bij_tear.';
                                Oa231=-H31.'+X1-Mgse_tear.'*Bij_tear.';
                                Oa241=-Dij_tear;
                                Oa251=lambda*H21;
                                Oa261=zeros(6,1);
                                Oa271=-J_tear*X1;
                                Oa281=kappa* Mgse_tear.'*Rbij.';
                                
                                Oa311=Oa131.';
                                Oa321=Oa231.';
                                Oa331=lambda*Q1+X1.'+X1;
                                Oa341=-Dij_tear;
                                Oa351=lambda*H31;
                                Oa361=zeros(6,1);
                                Oa371=-J_tear*X1;
                                Oa381=zeros(6,6);
                                
                                Oa411=Oa141.';
                                Oa421=Oa241.';
                                Oa431=Oa341.';
                                Oa441=-gamma*gamma*eye(1);
                                Oa451=zeros(1,6);
                                Oa461=zeros(1,1);
                                Oa471=zeros(1,6);
                                Oa481=zeros(1,6);
                                
                                Oa511=Oa151.';
                                Oa521=Oa251.';
                                Oa531=Oa351.';
                                Oa541=Oa451.';
                                Oa551=-lambda*Q1;
                                Oa561=zeros(6,1);
                                Oa571=zeros(6,6);
                                Oa581=zeros(6,6);
                                
                                
                                
                                Oa611=Oa161.';
                                Oa621=Oa261.';
                                Oa631=Oa361.';
                                Oa641=Oa461.';
                                Oa651=Oa561.';
                                Oa661=-eye(1);
                                Oa671=zeros(1,6);
                                Oa681=zeros(1,6);
                                
                                Oa711=Oa171.';
                                Oa721=Oa271.';
                                Oa731=Oa371.';
                                Oa741=Oa471.';
                                Oa751=Oa571.';
                                Oa761=Oa671.';
                                Oa771=-kappa*eye(6);
                                Oa781=zeros(6,6);
                                
                                Oa811=Oa181.';
                                Oa821=Oa281.';
                                Oa831=Oa381.';
                                Oa841=Oa481.';
                                Oa851=Oa581.';
                                Oa861=Oa681.';
                                Oa871=Oa781.';
                                Oa881=-kappa*eye(6);
                                
                                O22221=[Oa111  Oa121  Oa131  Oa141  Oa151  Oa161  Oa171  Oa181;
                                    Oa211  Oa221  Oa231  Oa241  Oa251  Oa261  Oa271  Oa281;
                                    Oa311  Oa321  Oa331  Oa341  Oa351  Oa361  Oa371  Oa381;
                                    Oa411  Oa421  Oa431  Oa441  Oa451  Oa461  Oa471  Oa481;
                                    Oa511  Oa521  Oa531  Oa541  Oa551  Oa561  Oa571  Oa581;
                                    Oa611  Oa621  Oa631  Oa641  Oa651  Oa661  Oa671  Oa681;
                                    Oa711  Oa721  Oa731  Oa741  Oa751  Oa761  Oa771  Oa781;
                                    Oa811  Oa821  Oa831  Oa841  Oa851  Oa861  Oa871  Oa881];
                                
                                
                                SOS4=-s3.'*([O22221]+eps3.*eye(38))*s3;
                                % �ŧi�������x�}
                                
                                prog = sosineq(prog,SOS4,'sparsemultipartite',{ vars.' s3.'});
                            end
                            %ijjj-1
                            for subblock=3
                                %% í�w����3-53,i<j=g=s,(4*3)x(4*3):
                                %�X��1222�M2111
                                % Aij_tear=inv(L1)*T*(mAi+nAj)*L1;
                                % Bij_tear=inv(L1)*T*(mAi+nAj);
                                % Gij_tear=mCi+nCj;
                                % Dij_tear=inv(L1)*T*(mDi+nDj);
                                % Mgse_tear=mMge+nMse;
                                
                                %1222�}�l
                                Aij_tear=inv(L1)*T*(m*A1+n*A2)*L1;
                                Bij_tear=inv(L1)*T*(m*B1+n*B2);
                                Cij_tear=(m*C1+n*C2)*L1;
                                Dij_tear=inv(L1)*T*(m*D1+n*D2);
                                Mgse_tear=m*M21+n*M21;
                                Sum=Sum11;
                                
                                Oa111=diff(P1,x2).*(A1(3,:)*x_h)+H11+H11.'-Aij_tear*X1-X1.'*Aij_tear+Sum;
                                Oa121=-H11+H21.'-Bij_tear*Mgse_tear-X1.'*Aij_tear.';
                                Oa131=P1+H31.'+X1-X1.'*Aij_tear;
                                Oa141=-Dij_tear;
                                Oa151=lambda*H11;
                                Oa161=(Cij_tear*X1).';
                                Oa171=-J_tear*X1;
                                Oa181=kappa*X1.'*Raij.';
                                
                                Oa211=Oa121.';
                                Oa221=-H21-H21.'-Bij_tear*Mgse_tear-Mgse_tear.'*Bij_tear.';
                                Oa231=-H31.'+X1-Mgse_tear.'*Bij_tear.';
                                Oa241=-Dij_tear;
                                Oa251=lambda*H21;
                                Oa261=zeros(6,1);
                                Oa271=-J_tear*X1;
                                Oa281=kappa* Mgse_tear.'*Rbij.';
                                
                                Oa311=Oa131.';
                                Oa321=Oa231.';
                                Oa331=lambda*Q1+X1.'+X1;
                                Oa341=-Dij_tear;
                                Oa351=lambda*H31;
                                Oa361=zeros(6,1);
                                Oa371=-J_tear*X1;
                                Oa381=zeros(6,6);
                                
                                Oa411=Oa141.';
                                Oa421=Oa241.';
                                Oa431=Oa341.';
                                Oa441=-gamma*gamma*eye(1);
                                Oa451=zeros(1,6);
                                Oa461=zeros(1,1);
                                Oa471=zeros(1,6);
                                Oa481=zeros(1,6);
                                
                                Oa511=Oa151.';
                                Oa521=Oa251.';
                                Oa531=Oa351.';
                                Oa541=Oa451.';
                                Oa551=-lambda*Q1;
                                Oa561=zeros(6,1);
                                Oa571=zeros(6,6);
                                Oa581=zeros(6,6);
                                
                                
                                
                                Oa611=Oa161.';
                                Oa621=Oa261.';
                                Oa631=Oa361.';
                                Oa641=Oa461.';
                                Oa651=Oa561.';
                                Oa661=-eye(1);
                                Oa671=zeros(1,6);
                                Oa681=zeros(1,6);
                                
                                Oa711=Oa171.';
                                Oa721=Oa271.';
                                Oa731=Oa371.';
                                Oa741=Oa471.';
                                Oa751=Oa571.';
                                Oa761=Oa671.';
                                Oa771=-kappa*eye(6);
                                Oa781=zeros(6,6);
                                
                                Oa811=Oa181.';
                                Oa821=Oa281.';
                                Oa831=Oa381.';
                                Oa841=Oa481.';
                                Oa851=Oa581.';
                                Oa861=Oa681.';
                                Oa871=Oa781.';
                                Oa881=-kappa*eye(6);
                                
                                O12221=[Oa111  Oa121  Oa131  Oa141  Oa151  Oa161  Oa171  Oa181;
                                    Oa211  Oa221  Oa231  Oa241  Oa251  Oa261  Oa271  Oa281;
                                    Oa311  Oa321  Oa331  Oa341  Oa351  Oa361  Oa371  Oa381;
                                    Oa411  Oa421  Oa431  Oa441  Oa451  Oa461  Oa471  Oa481;
                                    Oa511  Oa521  Oa531  Oa541  Oa551  Oa561  Oa571  Oa581;
                                    Oa611  Oa621  Oa631  Oa641  Oa651  Oa661  Oa671  Oa681;
                                    Oa711  Oa721  Oa731  Oa741  Oa751  Oa761  Oa771  Oa781;
                                    Oa811  Oa821  Oa831  Oa841  Oa851  Oa861  Oa871  Oa881];
                                
                                
                                %1222����
                                %2111�}�l
                                Aij_tear=inv(L1)*T*(m*A2+n*A1)*L1;
                                Bij_tear=inv(L1)*T*(m*B2+n*B1);
                                Cij_tear=(m*C2+n*C1)*L1;
                                Dij_tear=inv(L1)*T*(m*D2+n*D1);
                                Mgse_tear=m*M11+n*M11;
                                Sum=Sum21;
                                
                                Ob111=diff(P1,x2).*(A1(3,:)*x_h)+H11+H11.'-Aij_tear*X1-X1.'*Aij_tear+Sum;
                                Ob121=-H11+H21.'-Bij_tear*Mgse_tear-X1.'*Aij_tear.';
                                Ob131=P1+H31.'+X1-X1.'*Aij_tear;
                                Ob141=-Dij_tear;
                                Ob151=lambda*H11;
                                Ob161=(Cij_tear*X1).';
                                Ob171=-J_tear*X1;
                                Ob181=kappa*X1.'*Raij.';
                                
                                Ob211=Ob121.';
                                Ob221=-H21-H21.'-Bij_tear*Mgse_tear-Mgse_tear.'*Bij_tear.';
                                Ob231=-H31.'+X1-Mgse_tear.'*Bij_tear.';
                                Ob241=-Dij_tear;
                                Ob251=lambda*H21;
                                Ob261=zeros(6,1);
                                Ob271=-J_tear*X1;
                                Ob281=kappa* Mgse_tear.'*Rbij.';
                                
                                Ob311=Ob131.';
                                Ob321=Ob231.';
                                Ob331=lambda*Q1+X1.'+X1;
                                Ob341=-Dij_tear;
                                Ob351=lambda*H31;
                                Ob361=zeros(6,1);
                                Ob371=-J_tear*X1;
                                Ob381=zeros(6,6);
                                
                                Ob411=Ob141.';
                                Ob421=Ob241.';
                                Ob431=Ob341.';
                                Ob441=-gamma*gamma*eye(1);
                                Ob451=zeros(1,6);
                                Ob461=zeros(1,1);
                                Ob471=zeros(1,6);
                                Ob481=zeros(1,6);
                                
                                Ob511=Ob151.';
                                Ob521=Ob251.';
                                Ob531=Ob351.';
                                Ob541=Ob451.';
                                Ob551=-lambda*Q1;
                                Ob561=zeros(6,1);
                                Ob571=zeros(6,6);
                                Ob581=zeros(6,6);
                                
                                
                                
                                Ob611=Ob161.';
                                Ob621=Ob261.';
                                Ob631=Ob361.';
                                Ob641=Ob461.';
                                Ob651=Ob561.';
                                Ob661=-eye(1);
                                Ob671=zeros(1,6);
                                Ob681=zeros(1,6);
                                
                                Ob711=Ob171.';
                                Ob721=Ob271.';
                                Ob731=Ob371.';
                                Ob741=Ob471.';
                                Ob751=Ob571.';
                                Ob761=Ob671.';
                                Ob771=-kappa*eye(6);
                                Ob781=zeros(6,6);
                                
                                Ob811=Ob181.';
                                Ob821=Ob281.';
                                Ob831=Ob381.';
                                Ob841=Ob481.';
                                Ob851=Ob581.';
                                Ob861=Ob681.';
                                Ob871=Ob781.';
                                Ob881=-kappa*eye(6);
                                
                                O21111=[Ob111  Ob121  Ob131  Ob141  Ob151  Ob161  Ob171  Ob181;
                                    Ob211  Ob221  Ob231  Ob241  Ob251  Ob261  Ob271  Ob281;
                                    Ob311  Ob321  Ob331  Ob341  Ob351  Ob361  Ob371  Ob381;
                                    Ob411  Ob421  Ob431  Ob441  Ob451  Ob461  Ob471  Ob481;
                                    Ob511  Ob521  Ob531  Ob541  Ob551  Ob561  Ob571  Ob581;
                                    Ob611  Ob621  Ob631  Ob641  Ob651  Ob661  Ob671  Ob681;
                                    Ob711  Ob721  Ob731  Ob741  Ob751  Ob761  Ob771  Ob781;
                                    Ob811  Ob821  Ob831  Ob841  Ob851  Ob861  Ob871  Ob881];
                                
                                
                                
                                %2111����
                                SOS5=s3.'*(-[O12221]-[O21111]-eps3.*eye(38))*s3;
                                % �ŧi�������x�}
                                
                                prog = sosineq(prog,SOS5,'sparsemultipartite',{ vars.' s3.'});
                            end
                            %ijgg1
                            for subblock=4
                                %% í�w����3-53,i<j=g=s,(4*3)x(4*3):
                                %�]��i���N�A�����
                                %�X��1122�M2211
                                %�X��2122�M1211
                                % Aij_tear=inv(L1)*T*(mAi+nAj)*L1;
                                % Bij_tear=inv(L1)*T*(mAi+nAj);
                                % Gij_tear=mCi+nCj;
                                % Dij_tear=inv(L1)*T*(mDi+nDj);
                                % Mgse_tear=mMge+nMse;
                                %1122�}�l
                                Aij_tear=inv(L1)*T*(m*A1+n*A1)*L1;
                                Bij_tear=inv(L1)*T*(m*B1+n*B1);
                                Cij_tear=(m*C1+n*C1)*L1;
                                Dij_tear=inv(L1)*T*(m*D1+n*D1);
                                Mgse_tear=m*M21+n*M21;
                                Sum=Sum11;
                                
                                Oa111=diff(P1,x2).*(A1(3,:)*x_h)+H11+H11.'-Aij_tear*X1-X1.'*Aij_tear+Sum;
                                Oa121=-H11+H21.'-Bij_tear*Mgse_tear-X1.'*Aij_tear.';
                                Oa131=P1+H31.'+X1-X1.'*Aij_tear;
                                Oa141=-Dij_tear;
                                Oa151=lambda*H11;
                                Oa161=(Cij_tear*X1).';
                                Oa171=-J_tear*X1;
                                Oa181=kappa*X1.'*Raij.';
                                
                                Oa211=Oa121.';
                                Oa221=-H21-H21.'-Bij_tear*Mgse_tear-Mgse_tear.'*Bij_tear.';
                                Oa231=-H31.'+X1-Mgse_tear.'*Bij_tear.';
                                Oa241=-Dij_tear;
                                Oa251=lambda*H21;
                                Oa261=zeros(6,1);
                                Oa271=-J_tear*X1;
                                Oa281=kappa* Mgse_tear.'*Rbij.';
                                
                                Oa311=Oa131.';
                                Oa321=Oa231.';
                                Oa331=lambda*Q1+X1.'+X1;
                                Oa341=-Dij_tear;
                                Oa351=lambda*H31;
                                Oa361=zeros(6,1);
                                Oa371=-J_tear*X1;
                                Oa381=zeros(6,6);
                                
                                Oa411=Oa141.';
                                Oa421=Oa241.';
                                Oa431=Oa341.';
                                Oa441=-gamma*gamma*eye(1);
                                Oa451=zeros(1,6);
                                Oa461=zeros(1,1);
                                Oa471=zeros(1,6);
                                Oa481=zeros(1,6);
                                
                                Oa511=Oa151.';
                                Oa521=Oa251.';
                                Oa531=Oa351.';
                                Oa541=Oa451.';
                                Oa551=-lambda*Q1;
                                Oa561=zeros(6,1);
                                Oa571=zeros(6,6);
                                Oa581=zeros(6,6);
                                
                                
                                
                                Oa611=Oa161.';
                                Oa621=Oa261.';
                                Oa631=Oa361.';
                                Oa641=Oa461.';
                                Oa651=Oa561.';
                                Oa661=-eye(1);
                                Oa671=zeros(1,6);
                                Oa681=zeros(1,6);
                                
                                Oa711=Oa171.';
                                Oa721=Oa271.';
                                Oa731=Oa371.';
                                Oa741=Oa471.';
                                Oa751=Oa571.';
                                Oa761=Oa671.';
                                Oa771=-kappa*eye(6);
                                Oa781=zeros(6,6);
                                
                                Oa811=Oa181.';
                                Oa821=Oa281.';
                                Oa831=Oa381.';
                                Oa841=Oa481.';
                                Oa851=Oa581.';
                                Oa861=Oa681.';
                                Oa871=Oa781.';
                                Oa881=-kappa*eye(6);
                                
                                O11221=[Oa111  Oa121  Oa131  Oa141  Oa151  Oa161  Oa171  Oa181;
                                    Oa211  Oa221  Oa231  Oa241  Oa251  Oa261  Oa271  Oa281;
                                    Oa311  Oa321  Oa331  Oa341  Oa351  Oa361  Oa371  Oa381;
                                    Oa411  Oa421  Oa431  Oa441  Oa451  Oa461  Oa471  Oa481;
                                    Oa511  Oa521  Oa531  Oa541  Oa551  Oa561  Oa571  Oa581;
                                    Oa611  Oa621  Oa631  Oa641  Oa651  Oa661  Oa671  Oa681;
                                    Oa711  Oa721  Oa731  Oa741  Oa751  Oa761  Oa771  Oa781;
                                    Oa811  Oa821  Oa831  Oa841  Oa851  Oa861  Oa871  Oa881];
                                
                                %1122����
                                %2211�}�l
                                Aij_tear=inv(L1)*T*(m*A2+n*A2)*L1;
                                Bij_tear=inv(L1)*T*(m*B2+n*B2);
                                Cij_tear=(m*C2+n*C2)*L1;
                                Dij_tear=inv(L1)*T*(m*D2+n*D2);
                                Mgse_tear=m*M11+n*M11;
                                Sum=Sum21;
                                
                                Ob111=diff(P1,x2).*(A1(3,:)*x_h)+H11+H11.'-Aij_tear*X1-X1.'*Aij_tear+Sum;
                                Ob121=-H11+H21.'-Bij_tear*Mgse_tear-X1.'*Aij_tear.';
                                Ob131=P1+H31.'+X1-X1.'*Aij_tear;
                                Ob141=-Dij_tear;
                                Ob151=lambda*H11;
                                Ob161=(Cij_tear*X1).';
                                Ob171=-J_tear*X1;
                                Ob181=kappa*X1.'*Raij.';
                                
                                Ob211=Ob121.';
                                Ob221=-H21-H21.'-Bij_tear*Mgse_tear-Mgse_tear.'*Bij_tear.';
                                Ob231=-H31.'+X1-Mgse_tear.'*Bij_tear.';
                                Ob241=-Dij_tear;
                                Ob251=lambda*H21;
                                Ob261=zeros(6,1);
                                Ob271=-J_tear*X1;
                                Ob281=kappa* Mgse_tear.'*Rbij.';
                                
                                Ob311=Ob131.';
                                Ob321=Ob231.';
                                Ob331=lambda*Q1+X1.'+X1;
                                Ob341=-Dij_tear;
                                Ob351=lambda*H31;
                                Ob361=zeros(6,1);
                                Ob371=-J_tear*X1;
                                Ob381=zeros(6,6);
                                
                                Ob411=Ob141.';
                                Ob421=Ob241.';
                                Ob431=Ob341.';
                                Ob441=-gamma*gamma*eye(1);
                                Ob451=zeros(1,6);
                                Ob461=zeros(1,1);
                                Ob471=zeros(1,6);
                                Ob481=zeros(1,6);
                                
                                Ob511=Ob151.';
                                Ob521=Ob251.';
                                Ob531=Ob351.';
                                Ob541=Ob451.';
                                Ob551=-lambda*Q1;
                                Ob561=zeros(6,1);
                                Ob571=zeros(6,6);
                                Ob581=zeros(6,6);
                                
                                
                                
                                Ob611=Ob161.';
                                Ob621=Ob261.';
                                Ob631=Ob361.';
                                Ob641=Ob461.';
                                Ob651=Ob561.';
                                Ob661=-eye(1);
                                Ob671=zeros(1,6);
                                Ob681=zeros(1,6);
                                
                                Ob711=Ob171.';
                                Ob721=Ob271.';
                                Ob731=Ob371.';
                                Ob741=Ob471.';
                                Ob751=Ob571.';
                                Ob761=Ob671.';
                                Ob771=-kappa*eye(6);
                                Ob781=zeros(6,6);
                                
                                Ob811=Ob181.';
                                Ob821=Ob281.';
                                Ob831=Ob381.';
                                Ob841=Ob481.';
                                Ob851=Ob581.';
                                Ob861=Ob681.';
                                Ob871=Ob781.';
                                Ob881=-kappa*eye(6);
                                
                                O22111=[Ob111  Ob121  Ob131  Ob141  Ob151  Ob161  Ob171  Ob181;
                                    Ob211  Ob221  Ob231  Ob241  Ob251  Ob261  Ob271  Ob281;
                                    Ob311  Ob321  Ob331  Ob341  Ob351  Ob361  Ob371  Ob381;
                                    Ob411  Ob421  Ob431  Ob441  Ob451  Ob461  Ob471  Ob481;
                                    Ob511  Ob521  Ob531  Ob541  Ob551  Ob561  Ob571  Ob581;
                                    Ob611  Ob621  Ob631  Ob641  Ob651  Ob661  Ob671  Ob681;
                                    Ob711  Ob721  Ob731  Ob741  Ob751  Ob761  Ob771  Ob781;
                                    Ob811  Ob821  Ob831  Ob841  Ob851  Ob861  Ob871  Ob881];
                                
                                %2211����
                                SOS6=s3.'*(-[O11221]-[O22111]-eps3.*eye(38))*s3;
                                % �ŧi�������x�}
                                
                                prog = sosineq(prog,SOS6,'sparsemultipartite',{ vars.' s3.'});
                            end
                            for subblock=5
                                %% í�w����3-53,i<j=g=s,(4*3)x(4*3):
                                %�X��2122�M1211
                                % Aij_tear=inv(L1)*T*(mAi+nAj)*L1;
                                % Bij_tear=inv(L1)*T*(mAi+nAj);
                                % Gij_tear=mCi+nCj;
                                % Dij_tear=inv(L1)*T*(mDi+nDj);
                                % Mgse_tear=mMge+nMse;
                                %2122�}�l
                                Aij_tear=inv(L1)*T*(m*A2+n*A1)*L1;
                                Bij_tear=inv(L1)*T*(m*B2+n*B1);
                                Cij_tear=(m*C2+n*C1)*L1;
                                Dij_tear=inv(L1)*T*(m*D2+n*D1);
                                Mgse_tear=m*M21+n*M21;
                                Sum=Sum21;
                                
                                Oa111=diff(P1,x2).*(A1(3,:)*x_h)+H11+H11.'-Aij_tear*X1-X1.'*Aij_tear+Sum;
                                Oa121=-H11+H21.'-Bij_tear*Mgse_tear-X1.'*Aij_tear.';
                                Oa131=P1+H31.'+X1-X1.'*Aij_tear;
                                Oa141=-Dij_tear;
                                Oa151=lambda*H11;
                                Oa161=(Cij_tear*X1).';
                                Oa171=-J_tear*X1;
                                Oa181=kappa*X1.'*Raij.';
                                
                                Oa211=Oa121.';
                                Oa221=-H21-H21.'-Bij_tear*Mgse_tear-Mgse_tear.'*Bij_tear.';
                                Oa231=-H31.'+X1-Mgse_tear.'*Bij_tear.';
                                Oa241=-Dij_tear;
                                Oa251=lambda*H21;
                                Oa261=zeros(6,1);
                                Oa271=-J_tear*X1;
                                Oa281=kappa* Mgse_tear.'*Rbij.';
                                
                                Oa311=Oa131.';
                                Oa321=Oa231.';
                                Oa331=lambda*Q1+X1.'+X1;
                                Oa341=-Dij_tear;
                                Oa351=lambda*H31;
                                Oa361=zeros(6,1);
                                Oa371=-J_tear*X1;
                                Oa381=zeros(6,6);
                                
                                Oa411=Oa141.';
                                Oa421=Oa241.';
                                Oa431=Oa341.';
                                Oa441=-gamma*gamma*eye(1);
                                Oa451=zeros(1,6);
                                Oa461=zeros(1,1);
                                Oa471=zeros(1,6);
                                Oa481=zeros(1,6);
                                
                                Oa511=Oa151.';
                                Oa521=Oa251.';
                                Oa531=Oa351.';
                                Oa541=Oa451.';
                                Oa551=-lambda*Q1;
                                Oa561=zeros(6,1);
                                Oa571=zeros(6,6);
                                Oa581=zeros(6,6);
                                
                                
                                
                                Oa611=Oa161.';
                                Oa621=Oa261.';
                                Oa631=Oa361.';
                                Oa641=Oa461.';
                                Oa651=Oa561.';
                                Oa661=-eye(1);
                                Oa671=zeros(1,6);
                                Oa681=zeros(1,6);
                                
                                Oa711=Oa171.';
                                Oa721=Oa271.';
                                Oa731=Oa371.';
                                Oa741=Oa471.';
                                Oa751=Oa571.';
                                Oa761=Oa671.';
                                Oa771=-kappa*eye(6);
                                Oa781=zeros(6,6);
                                
                                Oa811=Oa181.';
                                Oa821=Oa281.';
                                Oa831=Oa381.';
                                Oa841=Oa481.';
                                Oa851=Oa581.';
                                Oa861=Oa681.';
                                Oa871=Oa781.';
                                Oa881=-kappa*eye(6);
                                
                                O21221=[Oa111  Oa121  Oa131  Oa141  Oa151  Oa161  Oa171  Oa181;
                                    Oa211  Oa221  Oa231  Oa241  Oa251  Oa261  Oa271  Oa281;
                                    Oa311  Oa321  Oa331  Oa341  Oa351  Oa361  Oa371  Oa381;
                                    Oa411  Oa421  Oa431  Oa441  Oa451  Oa461  Oa471  Oa481;
                                    Oa511  Oa521  Oa531  Oa541  Oa551  Oa561  Oa571  Oa581;
                                    Oa611  Oa621  Oa631  Oa641  Oa651  Oa661  Oa671  Oa681;
                                    Oa711  Oa721  Oa731  Oa741  Oa751  Oa761  Oa771  Oa781;
                                    Oa811  Oa821  Oa831  Oa841  Oa851  Oa861  Oa871  Oa881];
                                
                                %2122����
                                %1211�}�l
                                Aij_tear=inv(L1)*T*(m*A1+n*A2)*L1;
                                Bij_tear=inv(L1)*T*(m*B1+n*B2);
                                Cij_tear=(m*C1+n*C2)*L1;
                                Dij_tear=inv(L1)*T*(m*D1+n*D2);
                                Mgse_tear=m*M11+n*M11;
                                Sum=Sum11;
                                
                                Ob111=diff(P1,x2).*(A1(3,:)*x_h)+H11+H11.'-Aij_tear*X1-X1.'*Aij_tear+Sum;
                                Ob121=-H11+H21.'-Bij_tear*Mgse_tear-X1.'*Aij_tear.';
                                Ob131=P1+H31.'+X1-X1.'*Aij_tear;
                                Ob141=-Dij_tear;
                                Ob151=lambda*H11;
                                Ob161=(Cij_tear*X1).';
                                Ob171=-J_tear*X1;
                                Ob181=kappa*X1.'*Raij.';
                                
                                Ob211=Ob121.';
                                Ob221=-H21-H21.'-Bij_tear*Mgse_tear-Mgse_tear.'*Bij_tear.';
                                Ob231=-H31.'+X1-Mgse_tear.'*Bij_tear.';
                                Ob241=-Dij_tear;
                                Ob251=lambda*H21;
                                Ob261=zeros(6,1);
                                Ob271=-J_tear*X1;
                                Ob281=kappa* Mgse_tear.'*Rbij.';
                                
                                Ob311=Ob131.';
                                Ob321=Ob231.';
                                Ob331=lambda*Q1+X1.'+X1;
                                Ob341=-Dij_tear;
                                Ob351=lambda*H31;
                                Ob361=zeros(6,1);
                                Ob371=-J_tear*X1;
                                Ob381=zeros(6,6);
                                
                                Ob411=Ob141.';
                                Ob421=Ob241.';
                                Ob431=Ob341.';
                                Ob441=-gamma*gamma*eye(1);
                                Ob451=zeros(1,6);
                                Ob461=zeros(1,1);
                                Ob471=zeros(1,6);
                                Ob481=zeros(1,6);
                                
                                Ob511=Ob151.';
                                Ob521=Ob251.';
                                Ob531=Ob351.';
                                Ob541=Ob451.';
                                Ob551=-lambda*Q1;
                                Ob561=zeros(6,1);
                                Ob571=zeros(6,6);
                                Ob581=zeros(6,6);
                                
                                
                                
                                Ob611=Ob161.';
                                Ob621=Ob261.';
                                Ob631=Ob361.';
                                Ob641=Ob461.';
                                Ob651=Ob561.';
                                Ob661=-eye(1);
                                Ob671=zeros(1,6);
                                Ob681=zeros(1,6);
                                
                                Ob711=Ob171.';
                                Ob721=Ob271.';
                                Ob731=Ob371.';
                                Ob741=Ob471.';
                                Ob751=Ob571.';
                                Ob761=Ob671.';
                                Ob771=-kappa*eye(6);
                                Ob781=zeros(6,6);
                                
                                Ob811=Ob181.';
                                Ob821=Ob281.';
                                Ob831=Ob381.';
                                Ob841=Ob481.';
                                Ob851=Ob581.';
                                Ob861=Ob681.';
                                Ob871=Ob781.';
                                Ob881=-kappa*eye(6);
                                
                                O12111=[Ob111  Ob121  Ob131  Ob141  Ob151  Ob161  Ob171  Ob181;
                                    Ob211  Ob221  Ob231  Ob241  Ob251  Ob261  Ob271  Ob281;
                                    Ob311  Ob321  Ob331  Ob341  Ob351  Ob361  Ob371  Ob381;
                                    Ob411  Ob421  Ob431  Ob441  Ob451  Ob461  Ob471  Ob481;
                                    Ob511  Ob521  Ob531  Ob541  Ob551  Ob561  Ob571  Ob581;
                                    Ob611  Ob621  Ob631  Ob641  Ob651  Ob661  Ob671  Ob681;
                                    Ob711  Ob721  Ob731  Ob741  Ob751  Ob761  Ob771  Ob781;
                                    Ob811  Ob821  Ob831  Ob841  Ob851  Ob861  Ob871  Ob881];
                                
                                
                                %1211����
                                SOS7=s3.'*(-[O21221]-[O12111]-eps3.*eye(38))*s3;
                                % �ŧi�������x�}
                                
                                prog = sosineq(prog,SOS7,'sparsemultipartite',{ vars.' s3.'});
                            end
                            %ijgs1
                            for subblock=6
                                %% í�w����3-53,i<j=g=s,(4*3)x(4*3):
                                %�]��i,j���N�A���|��
                                %�X��1112�M2221
                                %�X��1212�M2121
                                %�X��2112�M1221
                                %�X��2212�M1121
                                % Aij_tear=inv(L1)*T*(mAi+nAj)*L1;
                                % Bij_tear=inv(L1)*T*(mAi+nAj);
                                % Gij_tear=mCi+nCj;
                                % Dij_tear=inv(L1)*T*(mDi+nDj);
                                % Mgse_tear=mMge+nMse;
                                %1112�}�l
                                Aij_tear=inv(L1)*T*(m*A1+n*A1)*L1;
                                Bij_tear=inv(L1)*T*(m*B1+n*B1);
                                Cij_tear=(m*C1+n*C1)*L1;
                                Dij_tear=inv(L1)*T*(m*D1+n*D1);
                                Mgse_tear=m*M11+n*M21;
                                Sum=Sum11;
                                
                                Oa111=diff(P1,x2).*(A1(3,:)*x_h)+H11+H11.'-Aij_tear*X1-X1.'*Aij_tear+Sum;
                                Oa121=-H11+H21.'-Bij_tear*Mgse_tear-X1.'*Aij_tear.';
                                Oa131=P1+H31.'+X1-X1.'*Aij_tear;
                                Oa141=-Dij_tear;
                                Oa151=lambda*H11;
                                Oa161=(Cij_tear*X1).';
                                Oa171=-J_tear*X1;
                                Oa181=kappa*X1.'*Raij.';
                                
                                Oa211=Oa121.';
                                Oa221=-H21-H21.'-Bij_tear*Mgse_tear-Mgse_tear.'*Bij_tear.';
                                Oa231=-H31.'+X1-Mgse_tear.'*Bij_tear.';
                                Oa241=-Dij_tear;
                                Oa251=lambda*H21;
                                Oa261=zeros(6,1);
                                Oa271=-J_tear*X1;
                                Oa281=kappa* Mgse_tear.'*Rbij.';
                                
                                Oa311=Oa131.';
                                Oa321=Oa231.';
                                Oa331=lambda*Q1+X1.'+X1;
                                Oa341=-Dij_tear;
                                Oa351=lambda*H31;
                                Oa361=zeros(6,1);
                                Oa371=-J_tear*X1;
                                Oa381=zeros(6,6);
                                
                                Oa411=Oa141.';
                                Oa421=Oa241.';
                                Oa431=Oa341.';
                                Oa441=-gamma*gamma*eye(1);
                                Oa451=zeros(1,6);
                                Oa461=zeros(1,1);
                                Oa471=zeros(1,6);
                                Oa481=zeros(1,6);
                                
                                Oa511=Oa151.';
                                Oa521=Oa251.';
                                Oa531=Oa351.';
                                Oa541=Oa451.';
                                Oa551=-lambda*Q1;
                                Oa561=zeros(6,1);
                                Oa571=zeros(6,6);
                                Oa581=zeros(6,6);
                                
                                
                                
                                Oa611=Oa161.';
                                Oa621=Oa261.';
                                Oa631=Oa361.';
                                Oa641=Oa461.';
                                Oa651=Oa561.';
                                Oa661=-eye(1);
                                Oa671=zeros(1,6);
                                Oa681=zeros(1,6);
                                
                                Oa711=Oa171.';
                                Oa721=Oa271.';
                                Oa731=Oa371.';
                                Oa741=Oa471.';
                                Oa751=Oa571.';
                                Oa761=Oa671.';
                                Oa771=-kappa*eye(6);
                                Oa781=zeros(6,6);
                                
                                Oa811=Oa181.';
                                Oa821=Oa281.';
                                Oa831=Oa381.';
                                Oa841=Oa481.';
                                Oa851=Oa581.';
                                Oa861=Oa681.';
                                Oa871=Oa781.';
                                Oa881=-kappa*eye(6);
                                
                                O11121=[Oa111  Oa121  Oa131  Oa141  Oa151  Oa161  Oa171  Oa181;
                                    Oa211  Oa221  Oa231  Oa241  Oa251  Oa261  Oa271  Oa281;
                                    Oa311  Oa321  Oa331  Oa341  Oa351  Oa361  Oa371  Oa381;
                                    Oa411  Oa421  Oa431  Oa441  Oa451  Oa461  Oa471  Oa481;
                                    Oa511  Oa521  Oa531  Oa541  Oa551  Oa561  Oa571  Oa581;
                                    Oa611  Oa621  Oa631  Oa641  Oa651  Oa661  Oa671  Oa681;
                                    Oa711  Oa721  Oa731  Oa741  Oa751  Oa761  Oa771  Oa781;
                                    Oa811  Oa821  Oa831  Oa841  Oa851  Oa861  Oa871  Oa881];
                                
                                
                                %1112����
                                %2221�}�l
                                Aij_tear=inv(L1)*T*(m*A2+n*A2)*L1;
                                Bij_tear=inv(L1)*T*(m*B2+n*B2);
                                Cij_tear=(m*C2+n*C2)*L1;
                                Dij_tear=inv(L1)*T*(m*D2+n*D2);
                                Mgse_tear=m*M21+n*M11;
                                Sum=Sum21;
                                
                                Ob111=diff(P1,x2).*(A1(3,:)*x_h)+H11+H11.'-Aij_tear*X1-X1.'*Aij_tear+Sum;
                                Ob121=-H11+H21.'-Bij_tear*Mgse_tear-X1.'*Aij_tear.';
                                Ob131=P1+H31.'+X1-X1.'*Aij_tear;
                                Ob141=-Dij_tear;
                                Ob151=lambda*H11;
                                Ob161=(Cij_tear*X1).';
                                Ob171=-J_tear*X1;
                                Ob181=kappa*X1.'*Raij.';
                                
                                Ob211=Ob121.';
                                Ob221=-H21-H21.'-Bij_tear*Mgse_tear-Mgse_tear.'*Bij_tear.';
                                Ob231=-H31.'+X1-Mgse_tear.'*Bij_tear.';
                                Ob241=-Dij_tear;
                                Ob251=lambda*H21;
                                Ob261=zeros(6,1);
                                Ob271=-J_tear*X1;
                                Ob281=kappa* Mgse_tear.'*Rbij.';
                                
                                Ob311=Ob131.';
                                Ob321=Ob231.';
                                Ob331=lambda*Q1+X1.'+X1;
                                Ob341=-Dij_tear;
                                Ob351=lambda*H31;
                                Ob361=zeros(6,1);
                                Ob371=-J_tear*X1;
                                Ob381=zeros(6,6);
                                
                                Ob411=Ob141.';
                                Ob421=Ob241.';
                                Ob431=Ob341.';
                                Ob441=-gamma*gamma*eye(1);
                                Ob451=zeros(1,6);
                                Ob461=zeros(1,1);
                                Ob471=zeros(1,6);
                                Ob481=zeros(1,6);
                                
                                Ob511=Ob151.';
                                Ob521=Ob251.';
                                Ob531=Ob351.';
                                Ob541=Ob451.';
                                Ob551=-lambda*Q1;
                                Ob561=zeros(6,1);
                                Ob571=zeros(6,6);
                                Ob581=zeros(6,6);
                                
                                
                                
                                Ob611=Ob161.';
                                Ob621=Ob261.';
                                Ob631=Ob361.';
                                Ob641=Ob461.';
                                Ob651=Ob561.';
                                Ob661=-eye(1);
                                Ob671=zeros(1,6);
                                Ob681=zeros(1,6);
                                
                                Ob711=Ob171.';
                                Ob721=Ob271.';
                                Ob731=Ob371.';
                                Ob741=Ob471.';
                                Ob751=Ob571.';
                                Ob761=Ob671.';
                                Ob771=-kappa*eye(6);
                                Ob781=zeros(6,6);
                                
                                Ob811=Ob181.';
                                Ob821=Ob281.';
                                Ob831=Ob381.';
                                Ob841=Ob481.';
                                Ob851=Ob581.';
                                Ob861=Ob681.';
                                Ob871=Ob781.';
                                Ob881=-kappa*eye(6);
                                
                                O22211=[Ob111  Ob121  Ob131  Ob141  Ob151  Ob161  Ob171  Ob181;
                                    Ob211  Ob221  Ob231  Ob241  Ob251  Ob261  Ob271  Ob281;
                                    Ob311  Ob321  Ob331  Ob341  Ob351  Ob361  Ob371  Ob381;
                                    Ob411  Ob421  Ob431  Ob441  Ob451  Ob461  Ob471  Ob481;
                                    Ob511  Ob521  Ob531  Ob541  Ob551  Ob561  Ob571  Ob581;
                                    Ob611  Ob621  Ob631  Ob641  Ob651  Ob661  Ob671  Ob681;
                                    Ob711  Ob721  Ob731  Ob741  Ob751  Ob761  Ob771  Ob781;
                                    Ob811  Ob821  Ob831  Ob841  Ob851  Ob861  Ob871  Ob881];
                                
                                %2221����
                                SOS8=s3.'*(-[O11121]-[O22211]-eps3.*eye(38))*s3;
                                % �ŧi�������x�}
                                
                                prog = sosineq(prog,SOS8,'sparsemultipartite',{ vars.' s3.'});
                            end
                            for subblock=7
                                %% í�w����3-53,i<j=g=s,(4*3)x(4*3):
                                %�X��1212�M2121
                                % Aij_tear=inv(L1)*T*(mAi+nAj)*L1;
                                % Bij_tear=inv(L1)*T*(mAi+nAj);
                                % Gij_tear=mCi+nCj;
                                % Dij_tear=inv(L1)*T*(mDi+nDj);
                                % Mgse_tear=mMge+nMse;
                                %1212�}�l
                                Aij_tear=inv(L1)*T*(m*A1+n*A2)*L1;
                                Bij_tear=inv(L1)*T*(m*B1+n*B2);
                                Cij_tear=(m*C1+n*C2)*L1;
                                Dij_tear=inv(L1)*T*(m*D1+n*D2);
                                Mgse_tear=m*M11+n*M21;
                                Sum=Sum11;
                                
                                Oa111=diff(P1,x2).*(A1(3,:)*x_h)+H11+H11.'-Aij_tear*X1-X1.'*Aij_tear+Sum;
                                Oa121=-H11+H21.'-Bij_tear*Mgse_tear-X1.'*Aij_tear.';
                                Oa131=P1+H31.'+X1-X1.'*Aij_tear;
                                Oa141=-Dij_tear;
                                Oa151=lambda*H11;
                                Oa161=(Cij_tear*X1).';
                                Oa171=-J_tear*X1;
                                Oa181=kappa*X1.'*Raij.';
                                
                                Oa211=Oa121.';
                                Oa221=-H21-H21.'-Bij_tear*Mgse_tear-Mgse_tear.'*Bij_tear.';
                                Oa231=-H31.'+X1-Mgse_tear.'*Bij_tear.';
                                Oa241=-Dij_tear;
                                Oa251=lambda*H21;
                                Oa261=zeros(6,1);
                                Oa271=-J_tear*X1;
                                Oa281=kappa* Mgse_tear.'*Rbij.';
                                
                                Oa311=Oa131.';
                                Oa321=Oa231.';
                                Oa331=lambda*Q1+X1.'+X1;
                                Oa341=-Dij_tear;
                                Oa351=lambda*H31;
                                Oa361=zeros(6,1);
                                Oa371=-J_tear*X1;
                                Oa381=zeros(6,6);
                                
                                Oa411=Oa141.';
                                Oa421=Oa241.';
                                Oa431=Oa341.';
                                Oa441=-gamma*gamma*eye(1);
                                Oa451=zeros(1,6);
                                Oa461=zeros(1,1);
                                Oa471=zeros(1,6);
                                Oa481=zeros(1,6);
                                
                                Oa511=Oa151.';
                                Oa521=Oa251.';
                                Oa531=Oa351.';
                                Oa541=Oa451.';
                                Oa551=-lambda*Q1;
                                Oa561=zeros(6,1);
                                Oa571=zeros(6,6);
                                Oa581=zeros(6,6);
                                
                                
                                
                                Oa611=Oa161.';
                                Oa621=Oa261.';
                                Oa631=Oa361.';
                                Oa641=Oa461.';
                                Oa651=Oa561.';
                                Oa661=-eye(1);
                                Oa671=zeros(1,6);
                                Oa681=zeros(1,6);
                                
                                Oa711=Oa171.';
                                Oa721=Oa271.';
                                Oa731=Oa371.';
                                Oa741=Oa471.';
                                Oa751=Oa571.';
                                Oa761=Oa671.';
                                Oa771=-kappa*eye(6);
                                Oa781=zeros(6,6);
                                
                                Oa811=Oa181.';
                                Oa821=Oa281.';
                                Oa831=Oa381.';
                                Oa841=Oa481.';
                                Oa851=Oa581.';
                                Oa861=Oa681.';
                                Oa871=Oa781.';
                                Oa881=-kappa*eye(6);
                                
                                O12121=[Oa111  Oa121  Oa131  Oa141  Oa151  Oa161  Oa171  Oa181;
                                    Oa211  Oa221  Oa231  Oa241  Oa251  Oa261  Oa271  Oa281;
                                    Oa311  Oa321  Oa331  Oa341  Oa351  Oa361  Oa371  Oa381;
                                    Oa411  Oa421  Oa431  Oa441  Oa451  Oa461  Oa471  Oa481;
                                    Oa511  Oa521  Oa531  Oa541  Oa551  Oa561  Oa571  Oa581;
                                    Oa611  Oa621  Oa631  Oa641  Oa651  Oa661  Oa671  Oa681;
                                    Oa711  Oa721  Oa731  Oa741  Oa751  Oa761  Oa771  Oa781;
                                    Oa811  Oa821  Oa831  Oa841  Oa851  Oa861  Oa871  Oa881];
                                
                                
                                %1212����
                                %2121�}�l
                                Aij_tear=inv(L1)*T*(m*A2+n*A1)*L1;
                                Bij_tear=inv(L1)*T*(m*B2+n*B1);
                                Cij_tear=(m*C2+n*C1)*L1;
                                Dij_tear=inv(L1)*T*(m*D2+n*D1);
                                Mgse_tear=m*M21+n*M11;
                                Sum=Sum21;
                                
                                Ob111=diff(P1,x2).*(A1(3,:)*x_h)+H11+H11.'-Aij_tear*X1-X1.'*Aij_tear+Sum;
                                Ob121=-H11+H21.'-Bij_tear*Mgse_tear-X1.'*Aij_tear.';
                                Ob131=P1+H31.'+X1-X1.'*Aij_tear;
                                Ob141=-Dij_tear;
                                Ob151=lambda*H11;
                                Ob161=(Cij_tear*X1).';
                                Ob171=-J_tear*X1;
                                Ob181=kappa*X1.'*Raij.';
                                
                                Ob211=Ob121.';
                                Ob221=-H21-H21.'-Bij_tear*Mgse_tear-Mgse_tear.'*Bij_tear.';
                                Ob231=-H31.'+X1-Mgse_tear.'*Bij_tear.';
                                Ob241=-Dij_tear;
                                Ob251=lambda*H21;
                                Ob261=zeros(6,1);
                                Ob271=-J_tear*X1;
                                Ob281=kappa* Mgse_tear.'*Rbij.';
                                
                                Ob311=Ob131.';
                                Ob321=Ob231.';
                                Ob331=lambda*Q1+X1.'+X1;
                                Ob341=-Dij_tear;
                                Ob351=lambda*H31;
                                Ob361=zeros(6,1);
                                Ob371=-J_tear*X1;
                                Ob381=zeros(6,6);
                                
                                Ob411=Ob141.';
                                Ob421=Ob241.';
                                Ob431=Ob341.';
                                Ob441=-gamma*gamma*eye(1);
                                Ob451=zeros(1,6);
                                Ob461=zeros(1,1);
                                Ob471=zeros(1,6);
                                Ob481=zeros(1,6);
                                
                                Ob511=Ob151.';
                                Ob521=Ob251.';
                                Ob531=Ob351.';
                                Ob541=Ob451.';
                                Ob551=-lambda*Q1;
                                Ob561=zeros(6,1);
                                Ob571=zeros(6,6);
                                Ob581=zeros(6,6);
                                
                                
                                
                                Ob611=Ob161.';
                                Ob621=Ob261.';
                                Ob631=Ob361.';
                                Ob641=Ob461.';
                                Ob651=Ob561.';
                                Ob661=-eye(1);
                                Ob671=zeros(1,6);
                                Ob681=zeros(1,6);
                                
                                Ob711=Ob171.';
                                Ob721=Ob271.';
                                Ob731=Ob371.';
                                Ob741=Ob471.';
                                Ob751=Ob571.';
                                Ob761=Ob671.';
                                Ob771=-kappa*eye(6);
                                Ob781=zeros(6,6);
                                
                                Ob811=Ob181.';
                                Ob821=Ob281.';
                                Ob831=Ob381.';
                                Ob841=Ob481.';
                                Ob851=Ob581.';
                                Ob861=Ob681.';
                                Ob871=Ob781.';
                                Ob881=-kappa*eye(6);
                                
                                O21211=[Ob111  Ob121  Ob131  Ob141  Ob151  Ob161  Ob171  Ob181;
                                    Ob211  Ob221  Ob231  Ob241  Ob251  Ob261  Ob271  Ob281;
                                    Ob311  Ob321  Ob331  Ob341  Ob351  Ob361  Ob371  Ob381;
                                    Ob411  Ob421  Ob431  Ob441  Ob451  Ob461  Ob471  Ob481;
                                    Ob511  Ob521  Ob531  Ob541  Ob551  Ob561  Ob571  Ob581;
                                    Ob611  Ob621  Ob631  Ob641  Ob651  Ob661  Ob671  Ob681;
                                    Ob711  Ob721  Ob731  Ob741  Ob751  Ob761  Ob771  Ob781;
                                    Ob811  Ob821  Ob831  Ob841  Ob851  Ob861  Ob871  Ob881];
                                
                                
                                %2121����
                                SOS9=s3.'*(-[O12121]-[O21211]-eps3.*eye(38))*s3;
                                % �ŧi�������x�}
                                
                                prog = sosineq(prog,SOS9,'sparsemultipartite',{ vars.' s3.'});
                            end
                            for subblock=8
                                %% í�w����3-53,i<j=g=s,(4*3)x(4*3):
                                %�X��2112�M1221
                                % Aij_tear=inv(L1)*T*(mAi+nAj)*L1;
                                % Bij_tear=inv(L1)*T*(mAi+nAj);
                                % Gij_tear=mCi+nCj;
                                % Dij_tear=inv(L1)*T*(mDi+nDj);
                                % Mgse_tear=mMge+nMse;
                                %2112�}�l
                                Aij_tear=inv(L1)*T*(m*A2+n*A1)*L1;
                                Bij_tear=inv(L1)*T*(m*B2+n*B1);
                                Cij_tear=(m*C2+n*C1)*L1;
                                Dij_tear=inv(L1)*T*(m*D2+n*D1);
                                Mgse_tear=m*M11+n*M21;
                                Sum=Sum21;
                                
                                Oa111=diff(P1,x2).*(A1(3,:)*x_h)+H11+H11.'-Aij_tear*X1-X1.'*Aij_tear+Sum;
                                Oa121=-H11+H21.'-Bij_tear*Mgse_tear-X1.'*Aij_tear.';
                                Oa131=P1+H31.'+X1-X1.'*Aij_tear;
                                Oa141=-Dij_tear;
                                Oa151=lambda*H11;
                                Oa161=(Cij_tear*X1).';
                                Oa171=-J_tear*X1;
                                Oa181=kappa*X1.'*Raij.';
                                
                                Oa211=Oa121.';
                                Oa221=-H21-H21.'-Bij_tear*Mgse_tear-Mgse_tear.'*Bij_tear.';
                                Oa231=-H31.'+X1-Mgse_tear.'*Bij_tear.';
                                Oa241=-Dij_tear;
                                Oa251=lambda*H21;
                                Oa261=zeros(6,1);
                                Oa271=-J_tear*X1;
                                Oa281=kappa* Mgse_tear.'*Rbij.';
                                
                                Oa311=Oa131.';
                                Oa321=Oa231.';
                                Oa331=lambda*Q1+X1.'+X1;
                                Oa341=-Dij_tear;
                                Oa351=lambda*H31;
                                Oa361=zeros(6,1);
                                Oa371=-J_tear*X1;
                                Oa381=zeros(6,6);
                                
                                Oa411=Oa141.';
                                Oa421=Oa241.';
                                Oa431=Oa341.';
                                Oa441=-gamma*gamma*eye(1);
                                Oa451=zeros(1,6);
                                Oa461=zeros(1,1);
                                Oa471=zeros(1,6);
                                Oa481=zeros(1,6);
                                
                                Oa511=Oa151.';
                                Oa521=Oa251.';
                                Oa531=Oa351.';
                                Oa541=Oa451.';
                                Oa551=-lambda*Q1;
                                Oa561=zeros(6,1);
                                Oa571=zeros(6,6);
                                Oa581=zeros(6,6);
                                
                                
                                
                                Oa611=Oa161.';
                                Oa621=Oa261.';
                                Oa631=Oa361.';
                                Oa641=Oa461.';
                                Oa651=Oa561.';
                                Oa661=-eye(1);
                                Oa671=zeros(1,6);
                                Oa681=zeros(1,6);
                                
                                Oa711=Oa171.';
                                Oa721=Oa271.';
                                Oa731=Oa371.';
                                Oa741=Oa471.';
                                Oa751=Oa571.';
                                Oa761=Oa671.';
                                Oa771=-kappa*eye(6);
                                Oa781=zeros(6,6);
                                
                                Oa811=Oa181.';
                                Oa821=Oa281.';
                                Oa831=Oa381.';
                                Oa841=Oa481.';
                                Oa851=Oa581.';
                                Oa861=Oa681.';
                                Oa871=Oa781.';
                                Oa881=-kappa*eye(6);
                                
                                O21121=[Oa111  Oa121  Oa131  Oa141  Oa151  Oa161  Oa171  Oa181;
                                    Oa211  Oa221  Oa231  Oa241  Oa251  Oa261  Oa271  Oa281;
                                    Oa311  Oa321  Oa331  Oa341  Oa351  Oa361  Oa371  Oa381;
                                    Oa411  Oa421  Oa431  Oa441  Oa451  Oa461  Oa471  Oa481;
                                    Oa511  Oa521  Oa531  Oa541  Oa551  Oa561  Oa571  Oa581;
                                    Oa611  Oa621  Oa631  Oa641  Oa651  Oa661  Oa671  Oa681;
                                    Oa711  Oa721  Oa731  Oa741  Oa751  Oa761  Oa771  Oa781;
                                    Oa811  Oa821  Oa831  Oa841  Oa851  Oa861  Oa871  Oa881];
                                
                                %2112����
                                %1221�}�l
                                Aij_tear=inv(L1)*T*(m*A1+n*A2)*L1;
                                Bij_tear=inv(L1)*T*(m*B1+n*B2);
                                Cij_tear=(m*C1+n*C2)*L1;
                                Dij_tear=inv(L1)*T*(m*D1+n*D2);
                                Mgse_tear=m*M21+n*M11;
                                Sum=Sum11;
                                
                                Ob111=diff(P1,x2).*(A1(3,:)*x_h)+H11+H11.'-Aij_tear*X1-X1.'*Aij_tear+Sum;
                                Ob121=-H11+H21.'-Bij_tear*Mgse_tear-X1.'*Aij_tear.';
                                Ob131=P1+H31.'+X1-X1.'*Aij_tear;
                                Ob141=-Dij_tear;
                                Ob151=lambda*H11;
                                Ob161=(Cij_tear*X1).';
                                Ob171=-J_tear*X1;
                                Ob181=kappa*X1.'*Raij.';
                                
                                Ob211=Ob121.';
                                Ob221=-H21-H21.'-Bij_tear*Mgse_tear-Mgse_tear.'*Bij_tear.';
                                Ob231=-H31.'+X1-Mgse_tear.'*Bij_tear.';
                                Ob241=-Dij_tear;
                                Ob251=lambda*H21;
                                Ob261=zeros(6,1);
                                Ob271=-J_tear*X1;
                                Ob281=kappa* Mgse_tear.'*Rbij.';
                                
                                Ob311=Ob131.';
                                Ob321=Ob231.';
                                Ob331=lambda*Q1+X1.'+X1;
                                Ob341=-Dij_tear;
                                Ob351=lambda*H31;
                                Ob361=zeros(6,1);
                                Ob371=-J_tear*X1;
                                Ob381=zeros(6,6);
                                
                                Ob411=Ob141.';
                                Ob421=Ob241.';
                                Ob431=Ob341.';
                                Ob441=-gamma*gamma*eye(1);
                                Ob451=zeros(1,6);
                                Ob461=zeros(1,1);
                                Ob471=zeros(1,6);
                                Ob481=zeros(1,6);
                                
                                Ob511=Ob151.';
                                Ob521=Ob251.';
                                Ob531=Ob351.';
                                Ob541=Ob451.';
                                Ob551=-lambda*Q1;
                                Ob561=zeros(6,1);
                                Ob571=zeros(6,6);
                                Ob581=zeros(6,6);
                                
                                
                                
                                Ob611=Ob161.';
                                Ob621=Ob261.';
                                Ob631=Ob361.';
                                Ob641=Ob461.';
                                Ob651=Ob561.';
                                Ob661=-eye(1);
                                Ob671=zeros(1,6);
                                Ob681=zeros(1,6);
                                
                                Ob711=Ob171.';
                                Ob721=Ob271.';
                                Ob731=Ob371.';
                                Ob741=Ob471.';
                                Ob751=Ob571.';
                                Ob761=Ob671.';
                                Ob771=-kappa*eye(6);
                                Ob781=zeros(6,6);
                                
                                Ob811=Ob181.';
                                Ob821=Ob281.';
                                Ob831=Ob381.';
                                Ob841=Ob481.';
                                Ob851=Ob581.';
                                Ob861=Ob681.';
                                Ob871=Ob781.';
                                Ob881=-kappa*eye(6);
                                
                                O12211=[Ob111  Ob121  Ob131  Ob141  Ob151  Ob161  Ob171  Ob181;
                                    Ob211  Ob221  Ob231  Ob241  Ob251  Ob261  Ob271  Ob281;
                                    Ob311  Ob321  Ob331  Ob341  Ob351  Ob361  Ob371  Ob381;
                                    Ob411  Ob421  Ob431  Ob441  Ob451  Ob461  Ob471  Ob481;
                                    Ob511  Ob521  Ob531  Ob541  Ob551  Ob561  Ob571  Ob581;
                                    Ob611  Ob621  Ob631  Ob641  Ob651  Ob661  Ob671  Ob681;
                                    Ob711  Ob721  Ob731  Ob741  Ob751  Ob761  Ob771  Ob781;
                                    Ob811  Ob821  Ob831  Ob841  Ob851  Ob861  Ob871  Ob881];
                                
                                
                                %1221����
                                SOS10=s3.'*(-[O21121]-[O12211]-eps3.*eye(38))*s3;
                                % �ŧi�������x�}
                                
                                prog = sosineq(prog,SOS10,'sparsemultipartite',{ vars.' s3.'});
                            end
                            for subblock=9
                                %% í�w����3-53,i<j=g=s,(4*3)x(4*3):
                                %�X��2212�M1121
                                % Aij_tear=inv(L1)*T*(mAi+nAj)*L1;
                                % Bij_tear=inv(L1)*T*(mAi+nAj);
                                % Gij_tear=mCi+nCj;
                                % Dij_tear=inv(L1)*T*(mDi+nDj);
                                % Mgse_tear=mMge+nMse;
                                %2212�}�l
                                Aij_tear=inv(L1)*T*(m*A2+n*A2)*L1;
                                Bij_tear=inv(L1)*T*(m*B2+n*B2);
                                Cij_tear=(m*C2+n*C2)*L1;
                                Dij_tear=inv(L1)*T*(m*D2+n*D2);
                                Mgse_tear=m*M11+n*M21;
                                Sum=Sum21;
                                
                                Oa111=diff(P1,x2).*(A1(3,:)*x_h)+H11+H11.'-Aij_tear*X1-X1.'*Aij_tear+Sum;
                                Oa121=-H11+H21.'-Bij_tear*Mgse_tear-X1.'*Aij_tear.';
                                Oa131=P1+H31.'+X1-X1.'*Aij_tear;
                                Oa141=-Dij_tear;
                                Oa151=lambda*H11;
                                Oa161=(Cij_tear*X1).';
                                Oa171=-J_tear*X1;
                                Oa181=kappa*X1.'*Raij.';
                                
                                Oa211=Oa121.';
                                Oa221=-H21-H21.'-Bij_tear*Mgse_tear-Mgse_tear.'*Bij_tear.';
                                Oa231=-H31.'+X1-Mgse_tear.'*Bij_tear.';
                                Oa241=-Dij_tear;
                                Oa251=lambda*H21;
                                Oa261=zeros(6,1);
                                Oa271=-J_tear*X1;
                                Oa281=kappa* Mgse_tear.'*Rbij.';
                                
                                Oa311=Oa131.';
                                Oa321=Oa231.';
                                Oa331=lambda*Q1+X1.'+X1;
                                Oa341=-Dij_tear;
                                Oa351=lambda*H31;
                                Oa361=zeros(6,1);
                                Oa371=-J_tear*X1;
                                Oa381=zeros(6,6);
                                
                                Oa411=Oa141.';
                                Oa421=Oa241.';
                                Oa431=Oa341.';
                                Oa441=-gamma*gamma*eye(1);
                                Oa451=zeros(1,6);
                                Oa461=zeros(1,1);
                                Oa471=zeros(1,6);
                                Oa481=zeros(1,6);
                                
                                Oa511=Oa151.';
                                Oa521=Oa251.';
                                Oa531=Oa351.';
                                Oa541=Oa451.';
                                Oa551=-lambda*Q1;
                                Oa561=zeros(6,1);
                                Oa571=zeros(6,6);
                                Oa581=zeros(6,6);
                                
                                
                                
                                Oa611=Oa161.';
                                Oa621=Oa261.';
                                Oa631=Oa361.';
                                Oa641=Oa461.';
                                Oa651=Oa561.';
                                Oa661=-eye(1);
                                Oa671=zeros(1,6);
                                Oa681=zeros(1,6);
                                
                                Oa711=Oa171.';
                                Oa721=Oa271.';
                                Oa731=Oa371.';
                                Oa741=Oa471.';
                                Oa751=Oa571.';
                                Oa761=Oa671.';
                                Oa771=-kappa*eye(6);
                                Oa781=zeros(6,6);
                                
                                Oa811=Oa181.';
                                Oa821=Oa281.';
                                Oa831=Oa381.';
                                Oa841=Oa481.';
                                Oa851=Oa581.';
                                Oa861=Oa681.';
                                Oa871=Oa781.';
                                Oa881=-kappa*eye(6);
                                
                                O22121=[Oa111  Oa121  Oa131  Oa141  Oa151  Oa161  Oa171  Oa181;
                                    Oa211  Oa221  Oa231  Oa241  Oa251  Oa261  Oa271  Oa281;
                                    Oa311  Oa321  Oa331  Oa341  Oa351  Oa361  Oa371  Oa381;
                                    Oa411  Oa421  Oa431  Oa441  Oa451  Oa461  Oa471  Oa481;
                                    Oa511  Oa521  Oa531  Oa541  Oa551  Oa561  Oa571  Oa581;
                                    Oa611  Oa621  Oa631  Oa641  Oa651  Oa661  Oa671  Oa681;
                                    Oa711  Oa721  Oa731  Oa741  Oa751  Oa761  Oa771  Oa781;
                                    Oa811  Oa821  Oa831  Oa841  Oa851  Oa861  Oa871  Oa881];
                                
                                %2212����
                                %1121�}�l
                                Aij_tear=inv(L1)*T*(m*A1+n*A1)*L1;
                                Bij_tear=inv(L1)*T*(m*B1+n*B1);
                                Cij_tear=(m*C1+n*C1)*L1;
                                Dij_tear=inv(L1)*T*(m*D1+n*D1);
                                Mgse_tear=m*M21+n*M11;
                                Sum=Sum11;
                                
                                Ob111=diff(P1,x2).*(A1(3,:)*x_h)+H11+H11.'-Aij_tear*X1-X1.'*Aij_tear+Sum;
                                Ob121=-H11+H21.'-Bij_tear*Mgse_tear-X1.'*Aij_tear.';
                                Ob131=P1+H31.'+X1-X1.'*Aij_tear;
                                Ob141=-Dij_tear;
                                Ob151=lambda*H11;
                                Ob161=(Cij_tear*X1).';
                                Ob171=-J_tear*X1;
                                Ob181=kappa*X1.'*Raij.';
                                
                                Ob211=Ob121.';
                                Ob221=-H21-H21.'-Bij_tear*Mgse_tear-Mgse_tear.'*Bij_tear.';
                                Ob231=-H31.'+X1-Mgse_tear.'*Bij_tear.';
                                Ob241=-Dij_tear;
                                Ob251=lambda*H21;
                                Ob261=zeros(6,1);
                                Ob271=-J_tear*X1;
                                Ob281=kappa* Mgse_tear.'*Rbij.';
                                
                                Ob311=Ob131.';
                                Ob321=Ob231.';
                                Ob331=lambda*Q1+X1.'+X1;
                                Ob341=-Dij_tear;
                                Ob351=lambda*H31;
                                Ob361=zeros(6,1);
                                Ob371=-J_tear*X1;
                                Ob381=zeros(6,6);
                                
                                Ob411=Ob141.';
                                Ob421=Ob241.';
                                Ob431=Ob341.';
                                Ob441=-gamma*gamma*eye(1);
                                Ob451=zeros(1,6);
                                Ob461=zeros(1,1);
                                Ob471=zeros(1,6);
                                Ob481=zeros(1,6);
                                
                                Ob511=Ob151.';
                                Ob521=Ob251.';
                                Ob531=Ob351.';
                                Ob541=Ob451.';
                                Ob551=-lambda*Q1;
                                Ob561=zeros(6,1);
                                Ob571=zeros(6,6);
                                Ob581=zeros(6,6);
                                
                                
                                
                                Ob611=Ob161.';
                                Ob621=Ob261.';
                                Ob631=Ob361.';
                                Ob641=Ob461.';
                                Ob651=Ob561.';
                                Ob661=-eye(1);
                                Ob671=zeros(1,6);
                                Ob681=zeros(1,6);
                                
                                Ob711=Ob171.';
                                Ob721=Ob271.';
                                Ob731=Ob371.';
                                Ob741=Ob471.';
                                Ob751=Ob571.';
                                Ob761=Ob671.';
                                Ob771=-kappa*eye(6);
                                Ob781=zeros(6,6);
                                
                                Ob811=Ob181.';
                                Ob821=Ob281.';
                                Ob831=Ob381.';
                                Ob841=Ob481.';
                                Ob851=Ob581.';
                                Ob861=Ob681.';
                                Ob871=Ob781.';
                                Ob881=-kappa*eye(6);
                                
                                O11211=[Ob111  Ob121  Ob131  Ob141  Ob151  Ob161  Ob171  Ob181;
                                    Ob211  Ob221  Ob231  Ob241  Ob251  Ob261  Ob271  Ob281;
                                    Ob311  Ob321  Ob331  Ob341  Ob351  Ob361  Ob371  Ob381;
                                    Ob411  Ob421  Ob431  Ob441  Ob451  Ob461  Ob471  Ob481;
                                    Ob511  Ob521  Ob531  Ob541  Ob551  Ob561  Ob571  Ob581;
                                    Ob611  Ob621  Ob631  Ob641  Ob651  Ob661  Ob671  Ob681;
                                    Ob711  Ob721  Ob731  Ob741  Ob751  Ob761  Ob771  Ob781;
                                    Ob811  Ob821  Ob831  Ob841  Ob851  Ob861  Ob871  Ob881];
                                
                                %2211����
                                SOS11=s3.'*(-[O22121]-[O11211]-eps3.*eye(38))*s3;
                                % �ŧi�������x�}
                                
                                prog = sosineq(prog,SOS11,'sparsemultipartite',{ vars.' s3.'});
                            end
                        end
                        for block=7
                            if N==2
                                Sum12 = 0;                                                              % �ŧiP1�x�}��ɶ��L��(i=1)���M
                                Sum22 = 0;                                                              % �ŧiP1�x�}��ɶ��L��(i=2)���M
                                
                                for i = 1:length(vars)                                                  % �̾�x_tilde���ƶq(i=1)�i��X��for�j��
                                    Sum12 = Sum12 + diff(P2,vars(i))*(TA1(i,:)*vars(1:sizeA));          % �p��P1�x�}��ɶ��L��(i=1)�ìۥ[
                                end
                                
                                for i = 1:length(vars)                                                  % �̾�x_tilde���ƶq(i=2)�i��X��for�j��
                                    Sum22 = Sum22 + diff(P2,vars(i))*(TA2(i,:)*vars(1:sizeA));          % �p��P1�x�}��ɶ��L��(i=2)�ìۥ[
                                end
                                
                                % �bSOS�{���[�J�x�}�����������A�������P2-eps1�x�}������
                                %eye(sizeA)�O�]��eps�O�`�ơA�x�}�����`��
                                %���w�z2í�w����Ĥ@��
                                SOS1 = s1.'*(P2-eps1.*eye(sizeA))*s1;
                                prog = sosmatrixineq(prog,SOS1);
                                % �bSOS�{���[�J�x�}�����������A�������Q2-eps2�x�}������
                                %���w�z2í�w����ĤG��
                                SOS2 = s2.'*(Q2-eps2.*eye(sizeA))*s2;
                                prog = sosmatrixineq(prog,SOS2);
                                
                                %N=2
                                
                                for subblock=1
                                    %% í�w����3-52,i=j=g=s=1,(4*3)x(4*3):
                                    %�]��i=j=g=s=1
                                    % Aij_tear=inv(L1)*T*(mAi+nAj)*L1;
                                    % Bij_tear=inv(L1)*T*(mAi+nAj);
                                    % Gij_tear=mCi+nCj;
                                    % Dij_tear=inv(L1)*T*(mDi+nDj);
                                    % Mgse_tear=mMge+nMse;
                                    Aij_tear=inv(L1)*T*(m*A1+n*A1)*L1;
                                    Bij_tear=inv(L1)*T*(m*B1+n*B1);
                                    Cij_tear=(m*C1+n*C1)*L1;
                                    Dij_tear=inv(L1)*T*(m*D1+n*D1);
                                    Mgse_tear=m*M12+n*M12;
                                    Sum=Sum12;
                                    
                                    Oa111=diff(P2,x2).*(A1(3,:)*x_h)+H12+H12.'-Aij_tear*X2-X2.'*Aij_tear+Sum;
                                    Oa121=-H12+H22.'-Bij_tear*Mgse_tear-X2.'*Aij_tear.';
                                    Oa131=P2+H32.'+X2-X2.'*Aij_tear;
                                    Oa141=-Dij_tear;
                                    Oa151=lambda*H12;
                                    Oa161=(Cij_tear*X2).';
                                    Oa171=-J_tear*X2;
                                    Oa181=kappa*X2.'*Raij.';
                                    
                                    Oa211=Oa121.';
                                    Oa221=-H22-H22.'-Bij_tear*Mgse_tear-Mgse_tear.'*Bij_tear.';
                                    Oa231=-H32.'+X2-Mgse_tear.'*Bij_tear.';
                                    Oa241=-Dij_tear;
                                    Oa251=lambda*H22;
                                    Oa261=zeros(6,1);
                                    Oa271=-J_tear*X2;
                                    Oa281=kappa* Mgse_tear.'*Rbij.';
                                    
                                    Oa311=Oa131.';
                                    Oa321=Oa231.';
                                    Oa331=lambda*Q2+X2.'+X2;
                                    Oa341=-Dij_tear;
                                    Oa351=lambda*H32;
                                    Oa361=zeros(6,1);
                                    Oa371=-J_tear*X2;
                                    Oa381=zeros(6,6);
                                    
                                    Oa411=Oa141.';
                                    Oa421=Oa241.';
                                    Oa431=Oa341.';
                                    Oa441=-gamma*gamma*eye(1);
                                    Oa451=zeros(1,6);
                                    Oa461=zeros(1,1);
                                    Oa471=zeros(1,6);
                                    Oa481=zeros(1,6);
                                    
                                    Oa511=Oa151.';
                                    Oa521=Oa251.';
                                    Oa531=Oa351.';
                                    Oa541=Oa451.';
                                    Oa551=-lambda*Q2;
                                    Oa561=zeros(6,1);
                                    Oa571=zeros(6,6);
                                    Oa581=zeros(6,6);
                                    
                                    
                                    
                                    Oa611=Oa161.';
                                    Oa621=Oa261.';
                                    Oa631=Oa361.';
                                    Oa641=Oa461.';
                                    Oa651=Oa561.';
                                    Oa661=-eye(1);
                                    Oa671=zeros(1,6);
                                    Oa681=zeros(1,6);
                                    
                                    Oa711=Oa171.';
                                    Oa721=Oa271.';
                                    Oa731=Oa371.';
                                    Oa741=Oa471.';
                                    Oa751=Oa571.';
                                    Oa761=Oa671.';
                                    Oa771=-kappa*eye(6);
                                    Oa781=zeros(6,6);
                                    
                                    Oa811=Oa181.';
                                    Oa821=Oa281.';
                                    Oa831=Oa381.';
                                    Oa841=Oa481.';
                                    Oa851=Oa581.';
                                    Oa861=Oa681.';
                                    Oa871=Oa781.';
                                    Oa881=-kappa*eye(6);
                                    
                                    
                                    
                                    O11111=[Oa111  Oa121  Oa131  Oa141  Oa151  Oa161  Oa171  Oa181;
                                        Oa211  Oa221  Oa231  Oa241  Oa251  Oa261  Oa271  Oa281;
                                        Oa311  Oa321  Oa331  Oa341  Oa351  Oa361  Oa371  Oa381;
                                        Oa411  Oa421  Oa431  Oa441  Oa451  Oa461  Oa471  Oa481;
                                        Oa511  Oa521  Oa531  Oa541  Oa551  Oa561  Oa571  Oa581;
                                        Oa611  Oa621  Oa631  Oa641  Oa651  Oa661  Oa671  Oa681;
                                        Oa711  Oa721  Oa731  Oa741  Oa751  Oa761  Oa771  Oa781;
                                        Oa811  Oa821  Oa831  Oa841  Oa851  Oa861  Oa871  Oa881];
                                    
                                    SOS3=-s3.'*([O11111]+eps3.*eye(38))*s3;
                                    % �ŧi�������x�}
                                    
                                    prog = sosineq(prog,SOS3,'sparsemultipartite',{ vars.' s3.'});
                                end
                                %22221
                                for subblock=2
                                    %% í�w����3-52,i=j=g=s=2,(4*3)x(4*3):
                                    %�]��i=j=g=s=2
                                    % Aij_tear=inv(L1)*T*(mAi+nAj)*L1;
                                    % Bij_tear=inv(L1)*T*(mAi+nAj);
                                    % Gij_tear=mCi+nCj;
                                    % Dij_tear=inv(L1)*T*(mDi+nDj);
                                    % Mgse_tear=mMge+nMse;
                                    Aij_tear=inv(L1)*T*(m*A2+n*A2)*L1;
                                    Bij_tear=inv(L1)*T*(m*B2+n*B2);
                                    Cij_tear=(m*C2+n*C2)*L1;
                                    Dij_tear=inv(L1)*T*(m*D2+n*D2);
                                    Mgse_tear=m*M22+n*M22;
                                    Sum=Sum22;
                                    
                                    Oa111=diff(P2,x2).*(A1(3,:)*x_h)+H12+H12.'-Aij_tear*X2-X2.'*Aij_tear+Sum;
                                    Oa121=-H12+H22.'-Bij_tear*Mgse_tear-X2.'*Aij_tear.';
                                    Oa131=P2+H32.'+X2-X2.'*Aij_tear;
                                    Oa141=-Dij_tear;
                                    Oa151=lambda*H12;
                                    Oa161=(Cij_tear*X2).';
                                    Oa171=-J_tear*X2;
                                    Oa181=kappa*X2.'*Raij.';
                                    
                                    Oa211=Oa121.';
                                    Oa221=-H22-H22.'-Bij_tear*Mgse_tear-Mgse_tear.'*Bij_tear.';
                                    Oa231=-H32.'+X2-Mgse_tear.'*Bij_tear.';
                                    Oa241=-Dij_tear;
                                    Oa251=lambda*H22;
                                    Oa261=zeros(6,1);
                                    Oa271=-J_tear*X2;
                                    Oa281=kappa* Mgse_tear.'*Rbij.';
                                    
                                    Oa311=Oa131.';
                                    Oa321=Oa231.';
                                    Oa331=lambda*Q2+X2.'+X2;
                                    Oa341=-Dij_tear;
                                    Oa351=lambda*H32;
                                    Oa361=zeros(6,1);
                                    Oa371=-J_tear*X2;
                                    Oa381=zeros(6,6);
                                    
                                    Oa411=Oa141.';
                                    Oa421=Oa241.';
                                    Oa431=Oa341.';
                                    Oa441=-gamma*gamma*eye(1);
                                    Oa451=zeros(1,6);
                                    Oa461=zeros(1,1);
                                    Oa471=zeros(1,6);
                                    Oa481=zeros(1,6);
                                    
                                    Oa511=Oa151.';
                                    Oa521=Oa251.';
                                    Oa531=Oa351.';
                                    Oa541=Oa451.';
                                    Oa551=-lambda*Q2;
                                    Oa561=zeros(6,1);
                                    Oa571=zeros(6,6);
                                    Oa581=zeros(6,6);
                                    
                                    
                                    
                                    Oa611=Oa161.';
                                    Oa621=Oa261.';
                                    Oa631=Oa361.';
                                    Oa641=Oa461.';
                                    Oa651=Oa561.';
                                    Oa661=-eye(1);
                                    Oa671=zeros(1,6);
                                    Oa681=zeros(1,6);
                                    
                                    Oa711=Oa171.';
                                    Oa721=Oa271.';
                                    Oa731=Oa371.';
                                    Oa741=Oa471.';
                                    Oa751=Oa571.';
                                    Oa761=Oa671.';
                                    Oa771=-kappa*eye(6);
                                    Oa781=zeros(6,6);
                                    
                                    Oa811=Oa181.';
                                    Oa821=Oa281.';
                                    Oa831=Oa381.';
                                    Oa841=Oa481.';
                                    Oa851=Oa581.';
                                    Oa861=Oa681.';
                                    Oa871=Oa781.';
                                    Oa881=-kappa*eye(6);
                                    
                                    O22221=[Oa111  Oa121  Oa131  Oa141  Oa151  Oa161  Oa171  Oa181;
                                        Oa211  Oa221  Oa231  Oa241  Oa251  Oa261  Oa271  Oa281;
                                        Oa311  Oa321  Oa331  Oa341  Oa351  Oa361  Oa371  Oa381;
                                        Oa411  Oa421  Oa431  Oa441  Oa451  Oa461  Oa471  Oa481;
                                        Oa511  Oa521  Oa531  Oa541  Oa551  Oa561  Oa571  Oa581;
                                        Oa611  Oa621  Oa631  Oa641  Oa651  Oa661  Oa671  Oa681;
                                        Oa711  Oa721  Oa731  Oa741  Oa751  Oa761  Oa771  Oa781;
                                        Oa811  Oa821  Oa831  Oa841  Oa851  Oa861  Oa871  Oa881];
                                    
                                    
                                    SOS4=-s3.'*([O22221]+eps3.*eye(38))*s3;
                                    % �ŧi�������x�}
                                    
                                    prog = sosineq(prog,SOS4,'sparsemultipartite',{ vars.' s3.'});
                                end
                                %ijjj-1
                                for subblock=3
                                    %% í�w����3-53,i<j=g=s,(4*3)x(4*3):
                                    %�X��1222�M2111
                                    % Aij_tear=inv(L1)*T*(mAi+nAj)*L1;
                                    % Bij_tear=inv(L1)*T*(mAi+nAj);
                                    % Gij_tear=mCi+nCj;
                                    % Dij_tear=inv(L1)*T*(mDi+nDj);
                                    % Mgse_tear=mMge+nMse;
                                    
                                    %1222�}�l
                                    Aij_tear=inv(L1)*T*(m*A1+n*A2)*L1;
                                    Bij_tear=inv(L1)*T*(m*B1+n*B2);
                                    Cij_tear=(m*C1+n*C2)*L1;
                                    Dij_tear=inv(L1)*T*(m*D1+n*D2);
                                    Mgse_tear=m*M22+n*M22;
                                    Sum=Sum12;
                                    
                                    Oa111=diff(P2,x2).*(A1(3,:)*x_h)+H12+H12.'-Aij_tear*X2-X2.'*Aij_tear+Sum;
                                    Oa121=-H12+H22.'-Bij_tear*Mgse_tear-X2.'*Aij_tear.';
                                    Oa131=P2+H32.'+X2-X2.'*Aij_tear;
                                    Oa141=-Dij_tear;
                                    Oa151=lambda*H12;
                                    Oa161=(Cij_tear*X2).';
                                    Oa171=-J_tear*X2;
                                    Oa181=kappa*X2.'*Raij.';
                                    
                                    Oa211=Oa121.';
                                    Oa221=-H22-H22.'-Bij_tear*Mgse_tear-Mgse_tear.'*Bij_tear.';
                                    Oa231=-H32.'+X2-Mgse_tear.'*Bij_tear.';
                                    Oa241=-Dij_tear;
                                    Oa251=lambda*H22;
                                    Oa261=zeros(6,1);
                                    Oa271=-J_tear*X2;
                                    Oa281=kappa* Mgse_tear.'*Rbij.';
                                    
                                    Oa311=Oa131.';
                                    Oa321=Oa231.';
                                    Oa331=lambda*Q2+X2.'+X2;
                                    Oa341=-Dij_tear;
                                    Oa351=lambda*H32;
                                    Oa361=zeros(6,1);
                                    Oa371=-J_tear*X2;
                                    Oa381=zeros(6,6);
                                    
                                    Oa411=Oa141.';
                                    Oa421=Oa241.';
                                    Oa431=Oa341.';
                                    Oa441=-gamma*gamma*eye(1);
                                    Oa451=zeros(1,6);
                                    Oa461=zeros(1,1);
                                    Oa471=zeros(1,6);
                                    Oa481=zeros(1,6);
                                    
                                    Oa511=Oa151.';
                                    Oa521=Oa251.';
                                    Oa531=Oa351.';
                                    Oa541=Oa451.';
                                    Oa551=-lambda*Q2;
                                    Oa561=zeros(6,1);
                                    Oa571=zeros(6,6);
                                    Oa581=zeros(6,6);
                                    
                                    
                                    
                                    Oa611=Oa161.';
                                    Oa621=Oa261.';
                                    Oa631=Oa361.';
                                    Oa641=Oa461.';
                                    Oa651=Oa561.';
                                    Oa661=-eye(1);
                                    Oa671=zeros(1,6);
                                    Oa681=zeros(1,6);
                                    
                                    Oa711=Oa171.';
                                    Oa721=Oa271.';
                                    Oa731=Oa371.';
                                    Oa741=Oa471.';
                                    Oa751=Oa571.';
                                    Oa761=Oa671.';
                                    Oa771=-kappa*eye(6);
                                    Oa781=zeros(6,6);
                                    
                                    Oa811=Oa181.';
                                    Oa821=Oa281.';
                                    Oa831=Oa381.';
                                    Oa841=Oa481.';
                                    Oa851=Oa581.';
                                    Oa861=Oa681.';
                                    Oa871=Oa781.';
                                    Oa881=-kappa*eye(6);
                                    
                                    O12221=[Oa111  Oa121  Oa131  Oa141  Oa151  Oa161  Oa171  Oa181;
                                        Oa211  Oa221  Oa231  Oa241  Oa251  Oa261  Oa271  Oa281;
                                        Oa311  Oa321  Oa331  Oa341  Oa351  Oa361  Oa371  Oa381;
                                        Oa411  Oa421  Oa431  Oa441  Oa451  Oa461  Oa471  Oa481;
                                        Oa511  Oa521  Oa531  Oa541  Oa551  Oa561  Oa571  Oa581;
                                        Oa611  Oa621  Oa631  Oa641  Oa651  Oa661  Oa671  Oa681;
                                        Oa711  Oa721  Oa731  Oa741  Oa751  Oa761  Oa771  Oa781;
                                        Oa811  Oa821  Oa831  Oa841  Oa851  Oa861  Oa871  Oa881];
                                    
                                    
                                    %1222����
                                    %2111�}�l
                                    Aij_tear=inv(L1)*T*(m*A2+n*A1)*L1;
                                    Bij_tear=inv(L1)*T*(m*B2+n*B1);
                                    Cij_tear=(m*C2+n*C1)*L1;
                                    Dij_tear=inv(L1)*T*(m*D2+n*D1);
                                    Mgse_tear=m*M12+n*M12;
                                    Sum=Sum22;
                                    
                                    Ob111=diff(P2,x2).*(A1(3,:)*x_h)+H12+H12.'-Aij_tear*X2-X2.'*Aij_tear+Sum;
                                    Ob121=-H12+H22.'-Bij_tear*Mgse_tear-X2.'*Aij_tear.';
                                    Ob131=P2+H32.'+X2-X2.'*Aij_tear;
                                    Ob141=-Dij_tear;
                                    Ob151=lambda*H12;
                                    Ob161=(Cij_tear*X2).';
                                    Ob171=-J_tear*X2;
                                    Ob181=kappa*X2.'*Raij.';
                                    
                                    Ob211=Ob121.';
                                    Ob221=-H22-H22.'-Bij_tear*Mgse_tear-Mgse_tear.'*Bij_tear.';
                                    Ob231=-H32.'+X2-Mgse_tear.'*Bij_tear.';
                                    Ob241=-Dij_tear;
                                    Ob251=lambda*H22;
                                    Ob261=zeros(6,1);
                                    Ob271=-J_tear*X2;
                                    Ob281=kappa* Mgse_tear.'*Rbij.';
                                    
                                    Ob311=Ob131.';
                                    Ob321=Ob231.';
                                    Ob331=lambda*Q2+X2.'+X2;
                                    Ob341=-Dij_tear;
                                    Ob351=lambda*H32;
                                    Ob361=zeros(6,1);
                                    Ob371=-J_tear*X2;
                                    Ob381=zeros(6,6);
                                    
                                    Ob411=Ob141.';
                                    Ob421=Ob241.';
                                    Ob431=Ob341.';
                                    Ob441=-gamma*gamma*eye(1);
                                    Ob451=zeros(1,6);
                                    Ob461=zeros(1,1);
                                    Ob471=zeros(1,6);
                                    Ob481=zeros(1,6);
                                    
                                    Ob511=Ob151.';
                                    Ob521=Ob251.';
                                    Ob531=Ob351.';
                                    Ob541=Ob451.';
                                    Ob551=-lambda*Q2;
                                    Ob561=zeros(6,1);
                                    Ob571=zeros(6,6);
                                    Ob581=zeros(6,6);
                                    
                                    
                                    
                                    Ob611=Ob161.';
                                    Ob621=Ob261.';
                                    Ob631=Ob361.';
                                    Ob641=Ob461.';
                                    Ob651=Ob561.';
                                    Ob661=-eye(1);
                                    Ob671=zeros(1,6);
                                    Ob681=zeros(1,6);
                                    
                                    Ob711=Ob171.';
                                    Ob721=Ob271.';
                                    Ob731=Ob371.';
                                    Ob741=Ob471.';
                                    Ob751=Ob571.';
                                    Ob761=Ob671.';
                                    Ob771=-kappa*eye(6);
                                    Ob781=zeros(6,6);
                                    
                                    Ob811=Ob181.';
                                    Ob821=Ob281.';
                                    Ob831=Ob381.';
                                    Ob841=Ob481.';
                                    Ob851=Ob581.';
                                    Ob861=Ob681.';
                                    Ob871=Ob781.';
                                    Ob881=-kappa*eye(6);
                                    
                                    O21111=[Ob111  Ob121  Ob131  Ob141  Ob151  Ob161  Ob171  Ob181;
                                        Ob211  Ob221  Ob231  Ob241  Ob251  Ob261  Ob271  Ob281;
                                        Ob311  Ob321  Ob331  Ob341  Ob351  Ob361  Ob371  Ob381;
                                        Ob411  Ob421  Ob431  Ob441  Ob451  Ob461  Ob471  Ob481;
                                        Ob511  Ob521  Ob531  Ob541  Ob551  Ob561  Ob571  Ob581;
                                        Ob611  Ob621  Ob631  Ob641  Ob651  Ob661  Ob671  Ob681;
                                        Ob711  Ob721  Ob731  Ob741  Ob751  Ob761  Ob771  Ob781;
                                        Ob811  Ob821  Ob831  Ob841  Ob851  Ob861  Ob871  Ob881];
                                    
                                    
                                    
                                    %2111����
                                    SOS5=s3.'*(-[O12221]-[O21111]-eps3.*eye(38))*s3;
                                    % �ŧi�������x�}
                                    
                                    prog = sosineq(prog,SOS5,'sparsemultipartite',{ vars.' s3.'});
                                end
                                %ijgg1
                                for subblock=4
                                    %% í�w����3-53,i<j=g=s,(4*3)x(4*3):
                                    %�]��i���N�A�����
                                    %�X��1122�M2211
                                    %�X��2122�M1211
                                    % Aij_tear=inv(L1)*T*(mAi+nAj)*L1;
                                    % Bij_tear=inv(L1)*T*(mAi+nAj);
                                    % Gij_tear=mCi+nCj;
                                    % Dij_tear=inv(L1)*T*(mDi+nDj);
                                    % Mgse_tear=mMge+nMse;
                                    %1122�}�l
                                    Aij_tear=inv(L1)*T*(m*A1+n*A1)*L1;
                                    Bij_tear=inv(L1)*T*(m*B1+n*B1);
                                    Cij_tear=(m*C1+n*C1)*L1;
                                    Dij_tear=inv(L1)*T*(m*D1+n*D1);
                                    Mgse_tear=m*M22+n*M22;
                                    Sum=Sum12;
                                    
                                    Oa111=diff(P2,x2).*(A1(3,:)*x_h)+H12+H12.'-Aij_tear*X2-X2.'*Aij_tear+Sum;
                                    Oa121=-H12+H22.'-Bij_tear*Mgse_tear-X2.'*Aij_tear.';
                                    Oa131=P2+H32.'+X2-X2.'*Aij_tear;
                                    Oa141=-Dij_tear;
                                    Oa151=lambda*H12;
                                    Oa161=(Cij_tear*X2).';
                                    Oa171=-J_tear*X2;
                                    Oa181=kappa*X2.'*Raij.';
                                    
                                    Oa211=Oa121.';
                                    Oa221=-H22-H22.'-Bij_tear*Mgse_tear-Mgse_tear.'*Bij_tear.';
                                    Oa231=-H32.'+X2-Mgse_tear.'*Bij_tear.';
                                    Oa241=-Dij_tear;
                                    Oa251=lambda*H22;
                                    Oa261=zeros(6,1);
                                    Oa271=-J_tear*X2;
                                    Oa281=kappa* Mgse_tear.'*Rbij.';
                                    
                                    Oa311=Oa131.';
                                    Oa321=Oa231.';
                                    Oa331=lambda*Q2+X2.'+X2;
                                    Oa341=-Dij_tear;
                                    Oa351=lambda*H32;
                                    Oa361=zeros(6,1);
                                    Oa371=-J_tear*X2;
                                    Oa381=zeros(6,6);
                                    
                                    Oa411=Oa141.';
                                    Oa421=Oa241.';
                                    Oa431=Oa341.';
                                    Oa441=-gamma*gamma*eye(1);
                                    Oa451=zeros(1,6);
                                    Oa461=zeros(1,1);
                                    Oa471=zeros(1,6);
                                    Oa481=zeros(1,6);
                                    
                                    Oa511=Oa151.';
                                    Oa521=Oa251.';
                                    Oa531=Oa351.';
                                    Oa541=Oa451.';
                                    Oa551=-lambda*Q2;
                                    Oa561=zeros(6,1);
                                    Oa571=zeros(6,6);
                                    Oa581=zeros(6,6);
                                    
                                    
                                    
                                    Oa611=Oa161.';
                                    Oa621=Oa261.';
                                    Oa631=Oa361.';
                                    Oa641=Oa461.';
                                    Oa651=Oa561.';
                                    Oa661=-eye(1);
                                    Oa671=zeros(1,6);
                                    Oa681=zeros(1,6);
                                    
                                    Oa711=Oa171.';
                                    Oa721=Oa271.';
                                    Oa731=Oa371.';
                                    Oa741=Oa471.';
                                    Oa751=Oa571.';
                                    Oa761=Oa671.';
                                    Oa771=-kappa*eye(6);
                                    Oa781=zeros(6,6);
                                    
                                    Oa811=Oa181.';
                                    Oa821=Oa281.';
                                    Oa831=Oa381.';
                                    Oa841=Oa481.';
                                    Oa851=Oa581.';
                                    Oa861=Oa681.';
                                    Oa871=Oa781.';
                                    Oa881=-kappa*eye(6);
                                    
                                    O11221=[Oa111  Oa121  Oa131  Oa141  Oa151  Oa161  Oa171  Oa181;
                                        Oa211  Oa221  Oa231  Oa241  Oa251  Oa261  Oa271  Oa281;
                                        Oa311  Oa321  Oa331  Oa341  Oa351  Oa361  Oa371  Oa381;
                                        Oa411  Oa421  Oa431  Oa441  Oa451  Oa461  Oa471  Oa481;
                                        Oa511  Oa521  Oa531  Oa541  Oa551  Oa561  Oa571  Oa581;
                                        Oa611  Oa621  Oa631  Oa641  Oa651  Oa661  Oa671  Oa681;
                                        Oa711  Oa721  Oa731  Oa741  Oa751  Oa761  Oa771  Oa781;
                                        Oa811  Oa821  Oa831  Oa841  Oa851  Oa861  Oa871  Oa881];
                                    
                                    %1122����
                                    %2211�}�l
                                    Aij_tear=inv(L1)*T*(m*A2+n*A2)*L1;
                                    Bij_tear=inv(L1)*T*(m*B2+n*B2);
                                    Cij_tear=(m*C2+n*C2)*L1;
                                    Dij_tear=inv(L1)*T*(m*D2+n*D2);
                                    Mgse_tear=m*M12+n*M12;
                                    Sum=Sum22;
                                    
                                    Ob111=diff(P2,x2).*(A1(3,:)*x_h)+H12+H12.'-Aij_tear*X2-X2.'*Aij_tear+Sum;
                                    Ob121=-H12+H22.'-Bij_tear*Mgse_tear-X2.'*Aij_tear.';
                                    Ob131=P2+H32.'+X2-X2.'*Aij_tear;
                                    Ob141=-Dij_tear;
                                    Ob151=lambda*H12;
                                    Ob161=(Cij_tear*X2).';
                                    Ob171=-J_tear*X2;
                                    Ob181=kappa*X2.'*Raij.';
                                    
                                    Ob211=Ob121.';
                                    Ob221=-H22-H22.'-Bij_tear*Mgse_tear-Mgse_tear.'*Bij_tear.';
                                    Ob231=-H32.'+X2-Mgse_tear.'*Bij_tear.';
                                    Ob241=-Dij_tear;
                                    Ob251=lambda*H22;
                                    Ob261=zeros(6,1);
                                    Ob271=-J_tear*X2;
                                    Ob281=kappa* Mgse_tear.'*Rbij.';
                                    
                                    Ob311=Ob131.';
                                    Ob321=Ob231.';
                                    Ob331=lambda*Q2+X2.'+X2;
                                    Ob341=-Dij_tear;
                                    Ob351=lambda*H32;
                                    Ob361=zeros(6,1);
                                    Ob371=-J_tear*X2;
                                    Ob381=zeros(6,6);
                                    
                                    Ob411=Ob141.';
                                    Ob421=Ob241.';
                                    Ob431=Ob341.';
                                    Ob441=-gamma*gamma*eye(1);
                                    Ob451=zeros(1,6);
                                    Ob461=zeros(1,1);
                                    Ob471=zeros(1,6);
                                    Ob481=zeros(1,6);
                                    
                                    Ob511=Ob151.';
                                    Ob521=Ob251.';
                                    Ob531=Ob351.';
                                    Ob541=Ob451.';
                                    Ob551=-lambda*Q2;
                                    Ob561=zeros(6,1);
                                    Ob571=zeros(6,6);
                                    Ob581=zeros(6,6);
                                    
                                    
                                    
                                    Ob611=Ob161.';
                                    Ob621=Ob261.';
                                    Ob631=Ob361.';
                                    Ob641=Ob461.';
                                    Ob651=Ob561.';
                                    Ob661=-eye(1);
                                    Ob671=zeros(1,6);
                                    Ob681=zeros(1,6);
                                    
                                    Ob711=Ob171.';
                                    Ob721=Ob271.';
                                    Ob731=Ob371.';
                                    Ob741=Ob471.';
                                    Ob751=Ob571.';
                                    Ob761=Ob671.';
                                    Ob771=-kappa*eye(6);
                                    Ob781=zeros(6,6);
                                    
                                    Ob811=Ob181.';
                                    Ob821=Ob281.';
                                    Ob831=Ob381.';
                                    Ob841=Ob481.';
                                    Ob851=Ob581.';
                                    Ob861=Ob681.';
                                    Ob871=Ob781.';
                                    Ob881=-kappa*eye(6);
                                    
                                    O22111=[Ob111  Ob121  Ob131  Ob141  Ob151  Ob161  Ob171  Ob181;
                                        Ob211  Ob221  Ob231  Ob241  Ob251  Ob261  Ob271  Ob281;
                                        Ob311  Ob321  Ob331  Ob341  Ob351  Ob361  Ob371  Ob381;
                                        Ob411  Ob421  Ob431  Ob441  Ob451  Ob461  Ob471  Ob481;
                                        Ob511  Ob521  Ob531  Ob541  Ob551  Ob561  Ob571  Ob581;
                                        Ob611  Ob621  Ob631  Ob641  Ob651  Ob661  Ob671  Ob681;
                                        Ob711  Ob721  Ob731  Ob741  Ob751  Ob761  Ob771  Ob781;
                                        Ob811  Ob821  Ob831  Ob841  Ob851  Ob861  Ob871  Ob881];
                                    
                                    %2211����
                                    SOS6=s3.'*(-[O11221]-[O22111]-eps3.*eye(38))*s3;
                                    % �ŧi�������x�}
                                    
                                    prog = sosineq(prog,SOS6,'sparsemultipartite',{ vars.' s3.'});
                                end
                                for subblock=5
                                    %% í�w����3-53,i<j=g=s,(4*3)x(4*3):
                                    %�X��2122�M1211
                                    % Aij_tear=inv(L1)*T*(mAi+nAj)*L1;
                                    % Bij_tear=inv(L1)*T*(mAi+nAj);
                                    % Gij_tear=mCi+nCj;
                                    % Dij_tear=inv(L1)*T*(mDi+nDj);
                                    % Mgse_tear=mMge+nMse;
                                    %2122�}�l
                                    Aij_tear=inv(L1)*T*(m*A2+n*A1)*L1;
                                    Bij_tear=inv(L1)*T*(m*B2+n*B1);
                                    Cij_tear=(m*C2+n*C1)*L1;
                                    Dij_tear=inv(L1)*T*(m*D2+n*D1);
                                    Mgse_tear=m*M22+n*M22;
                                    Sum=Sum22;
                                    
                                    Oa111=diff(P2,x2).*(A1(3,:)*x_h)+H12+H12.'-Aij_tear*X2-X2.'*Aij_tear+Sum;
                                    Oa121=-H12+H22.'-Bij_tear*Mgse_tear-X2.'*Aij_tear.';
                                    Oa131=P2+H32.'+X2-X2.'*Aij_tear;
                                    Oa141=-Dij_tear;
                                    Oa151=lambda*H12;
                                    Oa161=(Cij_tear*X2).';
                                    Oa171=-J_tear*X2;
                                    Oa181=kappa*X2.'*Raij.';
                                    
                                    Oa211=Oa121.';
                                    Oa221=-H22-H22.'-Bij_tear*Mgse_tear-Mgse_tear.'*Bij_tear.';
                                    Oa231=-H32.'+X2-Mgse_tear.'*Bij_tear.';
                                    Oa241=-Dij_tear;
                                    Oa251=lambda*H22;
                                    Oa261=zeros(6,1);
                                    Oa271=-J_tear*X2;
                                    Oa281=kappa* Mgse_tear.'*Rbij.';
                                    
                                    Oa311=Oa131.';
                                    Oa321=Oa231.';
                                    Oa331=lambda*Q2+X2.'+X2;
                                    Oa341=-Dij_tear;
                                    Oa351=lambda*H32;
                                    Oa361=zeros(6,1);
                                    Oa371=-J_tear*X2;
                                    Oa381=zeros(6,6);
                                    
                                    Oa411=Oa141.';
                                    Oa421=Oa241.';
                                    Oa431=Oa341.';
                                    Oa441=-gamma*gamma*eye(1);
                                    Oa451=zeros(1,6);
                                    Oa461=zeros(1,1);
                                    Oa471=zeros(1,6);
                                    Oa481=zeros(1,6);
                                    
                                    Oa511=Oa151.';
                                    Oa521=Oa251.';
                                    Oa531=Oa351.';
                                    Oa541=Oa451.';
                                    Oa551=-lambda*Q2;
                                    Oa561=zeros(6,1);
                                    Oa571=zeros(6,6);
                                    Oa581=zeros(6,6);
                                    
                                    
                                    
                                    Oa611=Oa161.';
                                    Oa621=Oa261.';
                                    Oa631=Oa361.';
                                    Oa641=Oa461.';
                                    Oa651=Oa561.';
                                    Oa661=-eye(1);
                                    Oa671=zeros(1,6);
                                    Oa681=zeros(1,6);
                                    
                                    Oa711=Oa171.';
                                    Oa721=Oa271.';
                                    Oa731=Oa371.';
                                    Oa741=Oa471.';
                                    Oa751=Oa571.';
                                    Oa761=Oa671.';
                                    Oa771=-kappa*eye(6);
                                    Oa781=zeros(6,6);
                                    
                                    Oa811=Oa181.';
                                    Oa821=Oa281.';
                                    Oa831=Oa381.';
                                    Oa841=Oa481.';
                                    Oa851=Oa581.';
                                    Oa861=Oa681.';
                                    Oa871=Oa781.';
                                    Oa881=-kappa*eye(6);
                                    
                                    O21221=[Oa111  Oa121  Oa131  Oa141  Oa151  Oa161  Oa171  Oa181;
                                        Oa211  Oa221  Oa231  Oa241  Oa251  Oa261  Oa271  Oa281;
                                        Oa311  Oa321  Oa331  Oa341  Oa351  Oa361  Oa371  Oa381;
                                        Oa411  Oa421  Oa431  Oa441  Oa451  Oa461  Oa471  Oa481;
                                        Oa511  Oa521  Oa531  Oa541  Oa551  Oa561  Oa571  Oa581;
                                        Oa611  Oa621  Oa631  Oa641  Oa651  Oa661  Oa671  Oa681;
                                        Oa711  Oa721  Oa731  Oa741  Oa751  Oa761  Oa771  Oa781;
                                        Oa811  Oa821  Oa831  Oa841  Oa851  Oa861  Oa871  Oa881];
                                    
                                    %2122����
                                    %1211�}�l
                                    Aij_tear=inv(L1)*T*(m*A1+n*A2)*L1;
                                    Bij_tear=inv(L1)*T*(m*B1+n*B2);
                                    Cij_tear=(m*C1+n*C2)*L1;
                                    Dij_tear=inv(L1)*T*(m*D1+n*D2);
                                    Mgse_tear=m*M12+n*M12;
                                    Sum=Sum12;
                                    
                                    Ob111=diff(P2,x2).*(A1(3,:)*x_h)+H12+H12.'-Aij_tear*X2-X2.'*Aij_tear+Sum;
                                    Ob121=-H12+H22.'-Bij_tear*Mgse_tear-X2.'*Aij_tear.';
                                    Ob131=P2+H32.'+X2-X2.'*Aij_tear;
                                    Ob141=-Dij_tear;
                                    Ob151=lambda*H12;
                                    Ob161=(Cij_tear*X2).';
                                    Ob171=-J_tear*X2;
                                    Ob181=kappa*X2.'*Raij.';
                                    
                                    Ob211=Ob121.';
                                    Ob221=-H22-H22.'-Bij_tear*Mgse_tear-Mgse_tear.'*Bij_tear.';
                                    Ob231=-H32.'+X2-Mgse_tear.'*Bij_tear.';
                                    Ob241=-Dij_tear;
                                    Ob251=lambda*H22;
                                    Ob261=zeros(6,1);
                                    Ob271=-J_tear*X2;
                                    Ob281=kappa* Mgse_tear.'*Rbij.';
                                    
                                    Ob311=Ob131.';
                                    Ob321=Ob231.';
                                    Ob331=lambda*Q2+X2.'+X2;
                                    Ob341=-Dij_tear;
                                    Ob351=lambda*H32;
                                    Ob361=zeros(6,1);
                                    Ob371=-J_tear*X2;
                                    Ob381=zeros(6,6);
                                    
                                    Ob411=Ob141.';
                                    Ob421=Ob241.';
                                    Ob431=Ob341.';
                                    Ob441=-gamma*gamma*eye(1);
                                    Ob451=zeros(1,6);
                                    Ob461=zeros(1,1);
                                    Ob471=zeros(1,6);
                                    Ob481=zeros(1,6);
                                    
                                    Ob511=Ob151.';
                                    Ob521=Ob251.';
                                    Ob531=Ob351.';
                                    Ob541=Ob451.';
                                    Ob551=-lambda*Q2;
                                    Ob561=zeros(6,1);
                                    Ob571=zeros(6,6);
                                    Ob581=zeros(6,6);
                                    
                                    
                                    
                                    Ob611=Ob161.';
                                    Ob621=Ob261.';
                                    Ob631=Ob361.';
                                    Ob641=Ob461.';
                                    Ob651=Ob561.';
                                    Ob661=-eye(1);
                                    Ob671=zeros(1,6);
                                    Ob681=zeros(1,6);
                                    
                                    Ob711=Ob171.';
                                    Ob721=Ob271.';
                                    Ob731=Ob371.';
                                    Ob741=Ob471.';
                                    Ob751=Ob571.';
                                    Ob761=Ob671.';
                                    Ob771=-kappa*eye(6);
                                    Ob781=zeros(6,6);
                                    
                                    Ob811=Ob181.';
                                    Ob821=Ob281.';
                                    Ob831=Ob381.';
                                    Ob841=Ob481.';
                                    Ob851=Ob581.';
                                    Ob861=Ob681.';
                                    Ob871=Ob781.';
                                    Ob881=-kappa*eye(6);
                                    
                                    O12111=[Ob111  Ob121  Ob131  Ob141  Ob151  Ob161  Ob171  Ob181;
                                        Ob211  Ob221  Ob231  Ob241  Ob251  Ob261  Ob271  Ob281;
                                        Ob311  Ob321  Ob331  Ob341  Ob351  Ob361  Ob371  Ob381;
                                        Ob411  Ob421  Ob431  Ob441  Ob451  Ob461  Ob471  Ob481;
                                        Ob511  Ob521  Ob531  Ob541  Ob551  Ob561  Ob571  Ob581;
                                        Ob611  Ob621  Ob631  Ob641  Ob651  Ob661  Ob671  Ob681;
                                        Ob711  Ob721  Ob731  Ob741  Ob751  Ob761  Ob771  Ob781;
                                        Ob811  Ob821  Ob831  Ob841  Ob851  Ob861  Ob871  Ob881];
                                    
                                    
                                    %1211����
                                    SOS7=s3.'*(-[O21221]-[O12111]-eps3.*eye(38))*s3;
                                    % �ŧi�������x�}
                                    
                                    prog = sosineq(prog,SOS7,'sparsemultipartite',{ vars.' s3.'});
                                end
                                %ijgs1
                                for subblock=6
                                    %% í�w����3-53,i<j=g=s,(4*3)x(4*3):
                                    %�]��i,j���N�A���|��
                                    %�X��1112�M2221
                                    %�X��1212�M2121
                                    %�X��2112�M1221
                                    %�X��2212�M1121
                                    % Aij_tear=inv(L1)*T*(mAi+nAj)*L1;
                                    % Bij_tear=inv(L1)*T*(mAi+nAj);
                                    % Gij_tear=mCi+nCj;
                                    % Dij_tear=inv(L1)*T*(mDi+nDj);
                                    % Mgse_tear=mMge+nMse;
                                    %1112�}�l
                                    Aij_tear=inv(L1)*T*(m*A1+n*A1)*L1;
                                    Bij_tear=inv(L1)*T*(m*B1+n*B1);
                                    Cij_tear=(m*C1+n*C1)*L1;
                                    Dij_tear=inv(L1)*T*(m*D1+n*D1);
                                    Mgse_tear=m*M12+n*M22;
                                    Sum=Sum12;
                                    
                                    Oa111=diff(P2,x2).*(A1(3,:)*x_h)+H12+H12.'-Aij_tear*X2-X2.'*Aij_tear+Sum;
                                    Oa121=-H12+H22.'-Bij_tear*Mgse_tear-X2.'*Aij_tear.';
                                    Oa131=P2+H32.'+X2-X2.'*Aij_tear;
                                    Oa141=-Dij_tear;
                                    Oa151=lambda*H12;
                                    Oa161=(Cij_tear*X2).';
                                    Oa171=-J_tear*X2;
                                    Oa181=kappa*X2.'*Raij.';
                                    
                                    Oa211=Oa121.';
                                    Oa221=-H22-H22.'-Bij_tear*Mgse_tear-Mgse_tear.'*Bij_tear.';
                                    Oa231=-H32.'+X2-Mgse_tear.'*Bij_tear.';
                                    Oa241=-Dij_tear;
                                    Oa251=lambda*H22;
                                    Oa261=zeros(6,1);
                                    Oa271=-J_tear*X2;
                                    Oa281=kappa* Mgse_tear.'*Rbij.';
                                    
                                    Oa311=Oa131.';
                                    Oa321=Oa231.';
                                    Oa331=lambda*Q2+X2.'+X2;
                                    Oa341=-Dij_tear;
                                    Oa351=lambda*H32;
                                    Oa361=zeros(6,1);
                                    Oa371=-J_tear*X2;
                                    Oa381=zeros(6,6);
                                    
                                    Oa411=Oa141.';
                                    Oa421=Oa241.';
                                    Oa431=Oa341.';
                                    Oa441=-gamma*gamma*eye(1);
                                    Oa451=zeros(1,6);
                                    Oa461=zeros(1,1);
                                    Oa471=zeros(1,6);
                                    Oa481=zeros(1,6);
                                    
                                    Oa511=Oa151.';
                                    Oa521=Oa251.';
                                    Oa531=Oa351.';
                                    Oa541=Oa451.';
                                    Oa551=-lambda*Q2;
                                    Oa561=zeros(6,1);
                                    Oa571=zeros(6,6);
                                    Oa581=zeros(6,6);
                                    
                                    
                                    
                                    Oa611=Oa161.';
                                    Oa621=Oa261.';
                                    Oa631=Oa361.';
                                    Oa641=Oa461.';
                                    Oa651=Oa561.';
                                    Oa661=-eye(1);
                                    Oa671=zeros(1,6);
                                    Oa681=zeros(1,6);
                                    
                                    Oa711=Oa171.';
                                    Oa721=Oa271.';
                                    Oa731=Oa371.';
                                    Oa741=Oa471.';
                                    Oa751=Oa571.';
                                    Oa761=Oa671.';
                                    Oa771=-kappa*eye(6);
                                    Oa781=zeros(6,6);
                                    
                                    Oa811=Oa181.';
                                    Oa821=Oa281.';
                                    Oa831=Oa381.';
                                    Oa841=Oa481.';
                                    Oa851=Oa581.';
                                    Oa861=Oa681.';
                                    Oa871=Oa781.';
                                    Oa881=-kappa*eye(6);
                                    
                                    O11121=[Oa111  Oa121  Oa131  Oa141  Oa151  Oa161  Oa171  Oa181;
                                        Oa211  Oa221  Oa231  Oa241  Oa251  Oa261  Oa271  Oa281;
                                        Oa311  Oa321  Oa331  Oa341  Oa351  Oa361  Oa371  Oa381;
                                        Oa411  Oa421  Oa431  Oa441  Oa451  Oa461  Oa471  Oa481;
                                        Oa511  Oa521  Oa531  Oa541  Oa551  Oa561  Oa571  Oa581;
                                        Oa611  Oa621  Oa631  Oa641  Oa651  Oa661  Oa671  Oa681;
                                        Oa711  Oa721  Oa731  Oa741  Oa751  Oa761  Oa771  Oa781;
                                        Oa811  Oa821  Oa831  Oa841  Oa851  Oa861  Oa871  Oa881];
                                    
                                    
                                    %1112����
                                    %2221�}�l
                                    Aij_tear=inv(L1)*T*(m*A2+n*A2)*L1;
                                    Bij_tear=inv(L1)*T*(m*B2+n*B2);
                                    Cij_tear=(m*C2+n*C2)*L1;
                                    Dij_tear=inv(L1)*T*(m*D2+n*D2);
                                    Mgse_tear=m*M22+n*M12;
                                    Sum=Sum22;
                                    
                                    Ob111=diff(P2,x2).*(A1(3,:)*x_h)+H12+H12.'-Aij_tear*X2-X2.'*Aij_tear+Sum;
                                    Ob121=-H12+H22.'-Bij_tear*Mgse_tear-X2.'*Aij_tear.';
                                    Ob131=P2+H32.'+X2-X2.'*Aij_tear;
                                    Ob141=-Dij_tear;
                                    Ob151=lambda*H12;
                                    Ob161=(Cij_tear*X2).';
                                    Ob171=-J_tear*X2;
                                    Ob181=kappa*X2.'*Raij.';
                                    
                                    Ob211=Ob121.';
                                    Ob221=-H22-H22.'-Bij_tear*Mgse_tear-Mgse_tear.'*Bij_tear.';
                                    Ob231=-H32.'+X2-Mgse_tear.'*Bij_tear.';
                                    Ob241=-Dij_tear;
                                    Ob251=lambda*H22;
                                    Ob261=zeros(6,1);
                                    Ob271=-J_tear*X2;
                                    Ob281=kappa* Mgse_tear.'*Rbij.';
                                    
                                    Ob311=Ob131.';
                                    Ob321=Ob231.';
                                    Ob331=lambda*Q2+X2.'+X2;
                                    Ob341=-Dij_tear;
                                    Ob351=lambda*H32;
                                    Ob361=zeros(6,1);
                                    Ob371=-J_tear*X2;
                                    Ob381=zeros(6,6);
                                    
                                    Ob411=Ob141.';
                                    Ob421=Ob241.';
                                    Ob431=Ob341.';
                                    Ob441=-gamma*gamma*eye(1);
                                    Ob451=zeros(1,6);
                                    Ob461=zeros(1,1);
                                    Ob471=zeros(1,6);
                                    Ob481=zeros(1,6);
                                    
                                    Ob511=Ob151.';
                                    Ob521=Ob251.';
                                    Ob531=Ob351.';
                                    Ob541=Ob451.';
                                    Ob551=-lambda*Q2;
                                    Ob561=zeros(6,1);
                                    Ob571=zeros(6,6);
                                    Ob581=zeros(6,6);
                                    
                                    
                                    
                                    Ob611=Ob161.';
                                    Ob621=Ob261.';
                                    Ob631=Ob361.';
                                    Ob641=Ob461.';
                                    Ob651=Ob561.';
                                    Ob661=-eye(1);
                                    Ob671=zeros(1,6);
                                    Ob681=zeros(1,6);
                                    
                                    Ob711=Ob171.';
                                    Ob721=Ob271.';
                                    Ob731=Ob371.';
                                    Ob741=Ob471.';
                                    Ob751=Ob571.';
                                    Ob761=Ob671.';
                                    Ob771=-kappa*eye(6);
                                    Ob781=zeros(6,6);
                                    
                                    Ob811=Ob181.';
                                    Ob821=Ob281.';
                                    Ob831=Ob381.';
                                    Ob841=Ob481.';
                                    Ob851=Ob581.';
                                    Ob861=Ob681.';
                                    Ob871=Ob781.';
                                    Ob881=-kappa*eye(6);
                                    
                                    O22211=[Ob111  Ob121  Ob131  Ob141  Ob151  Ob161  Ob171  Ob181;
                                        Ob211  Ob221  Ob231  Ob241  Ob251  Ob261  Ob271  Ob281;
                                        Ob311  Ob321  Ob331  Ob341  Ob351  Ob361  Ob371  Ob381;
                                        Ob411  Ob421  Ob431  Ob441  Ob451  Ob461  Ob471  Ob481;
                                        Ob511  Ob521  Ob531  Ob541  Ob551  Ob561  Ob571  Ob581;
                                        Ob611  Ob621  Ob631  Ob641  Ob651  Ob661  Ob671  Ob681;
                                        Ob711  Ob721  Ob731  Ob741  Ob751  Ob761  Ob771  Ob781;
                                        Ob811  Ob821  Ob831  Ob841  Ob851  Ob861  Ob871  Ob881];
                                    
                                    %2221����
                                    SOS8=s3.'*(-[O11121]-[O22211]-eps3.*eye(38))*s3;
                                    % �ŧi�������x�}
                                    
                                    prog = sosineq(prog,SOS8,'sparsemultipartite',{ vars.' s3.'});
                                end
                                for subblock=7
                                    %% í�w����3-53,i<j=g=s,(4*3)x(4*3):
                                    %�X��1212�M2121
                                    % Aij_tear=inv(L1)*T*(mAi+nAj)*L1;
                                    % Bij_tear=inv(L1)*T*(mAi+nAj);
                                    % Gij_tear=mCi+nCj;
                                    % Dij_tear=inv(L1)*T*(mDi+nDj);
                                    % Mgse_tear=mMge+nMse;
                                    %1212�}�l
                                    Aij_tear=inv(L1)*T*(m*A1+n*A2)*L1;
                                    Bij_tear=inv(L1)*T*(m*B1+n*B2);
                                    Cij_tear=(m*C1+n*C2)*L1;
                                    Dij_tear=inv(L1)*T*(m*D1+n*D2);
                                    Mgse_tear=m*M12+n*M22;
                                    Sum=Sum12;
                                    
                                    Oa111=diff(P2,x2).*(A1(3,:)*x_h)+H12+H12.'-Aij_tear*X2-X2.'*Aij_tear+Sum;
                                    Oa121=-H12+H22.'-Bij_tear*Mgse_tear-X2.'*Aij_tear.';
                                    Oa131=P2+H32.'+X2-X2.'*Aij_tear;
                                    Oa141=-Dij_tear;
                                    Oa151=lambda*H12;
                                    Oa161=(Cij_tear*X2).';
                                    Oa171=-J_tear*X2;
                                    Oa181=kappa*X2.'*Raij.';
                                    
                                    Oa211=Oa121.';
                                    Oa221=-H22-H22.'-Bij_tear*Mgse_tear-Mgse_tear.'*Bij_tear.';
                                    Oa231=-H32.'+X2-Mgse_tear.'*Bij_tear.';
                                    Oa241=-Dij_tear;
                                    Oa251=lambda*H22;
                                    Oa261=zeros(6,1);
                                    Oa271=-J_tear*X2;
                                    Oa281=kappa* Mgse_tear.'*Rbij.';
                                    
                                    Oa311=Oa131.';
                                    Oa321=Oa231.';
                                    Oa331=lambda*Q2+X2.'+X2;
                                    Oa341=-Dij_tear;
                                    Oa351=lambda*H32;
                                    Oa361=zeros(6,1);
                                    Oa371=-J_tear*X2;
                                    Oa381=zeros(6,6);
                                    
                                    Oa411=Oa141.';
                                    Oa421=Oa241.';
                                    Oa431=Oa341.';
                                    Oa441=-gamma*gamma*eye(1);
                                    Oa451=zeros(1,6);
                                    Oa461=zeros(1,1);
                                    Oa471=zeros(1,6);
                                    Oa481=zeros(1,6);
                                    
                                    Oa511=Oa151.';
                                    Oa521=Oa251.';
                                    Oa531=Oa351.';
                                    Oa541=Oa451.';
                                    Oa551=-lambda*Q2;
                                    Oa561=zeros(6,1);
                                    Oa571=zeros(6,6);
                                    Oa581=zeros(6,6);
                                    
                                    
                                    
                                    Oa611=Oa161.';
                                    Oa621=Oa261.';
                                    Oa631=Oa361.';
                                    Oa641=Oa461.';
                                    Oa651=Oa561.';
                                    Oa661=-eye(1);
                                    Oa671=zeros(1,6);
                                    Oa681=zeros(1,6);
                                    
                                    Oa711=Oa171.';
                                    Oa721=Oa271.';
                                    Oa731=Oa371.';
                                    Oa741=Oa471.';
                                    Oa751=Oa571.';
                                    Oa761=Oa671.';
                                    Oa771=-kappa*eye(6);
                                    Oa781=zeros(6,6);
                                    
                                    Oa811=Oa181.';
                                    Oa821=Oa281.';
                                    Oa831=Oa381.';
                                    Oa841=Oa481.';
                                    Oa851=Oa581.';
                                    Oa861=Oa681.';
                                    Oa871=Oa781.';
                                    Oa881=-kappa*eye(6);
                                    
                                    O12121=[Oa111  Oa121  Oa131  Oa141  Oa151  Oa161  Oa171  Oa181;
                                        Oa211  Oa221  Oa231  Oa241  Oa251  Oa261  Oa271  Oa281;
                                        Oa311  Oa321  Oa331  Oa341  Oa351  Oa361  Oa371  Oa381;
                                        Oa411  Oa421  Oa431  Oa441  Oa451  Oa461  Oa471  Oa481;
                                        Oa511  Oa521  Oa531  Oa541  Oa551  Oa561  Oa571  Oa581;
                                        Oa611  Oa621  Oa631  Oa641  Oa651  Oa661  Oa671  Oa681;
                                        Oa711  Oa721  Oa731  Oa741  Oa751  Oa761  Oa771  Oa781;
                                        Oa811  Oa821  Oa831  Oa841  Oa851  Oa861  Oa871  Oa881];
                                    
                                    
                                    %1212����
                                    %2121�}�l
                                    Aij_tear=inv(L1)*T*(m*A2+n*A1)*L1;
                                    Bij_tear=inv(L1)*T*(m*B2+n*B1);
                                    Cij_tear=(m*C2+n*C1)*L1;
                                    Dij_tear=inv(L1)*T*(m*D2+n*D1);
                                    Mgse_tear=m*M22+n*M12;
                                    Sum=Sum22;
                                    
                                    Ob111=diff(P2,x2).*(A1(3,:)*x_h)+H12+H12.'-Aij_tear*X2-X2.'*Aij_tear+Sum;
                                    Ob121=-H12+H22.'-Bij_tear*Mgse_tear-X2.'*Aij_tear.';
                                    Ob131=P2+H32.'+X2-X2.'*Aij_tear;
                                    Ob141=-Dij_tear;
                                    Ob151=lambda*H12;
                                    Ob161=(Cij_tear*X2).';
                                    Ob171=-J_tear*X2;
                                    Ob181=kappa*X2.'*Raij.';
                                    
                                    Ob211=Ob121.';
                                    Ob221=-H22-H22.'-Bij_tear*Mgse_tear-Mgse_tear.'*Bij_tear.';
                                    Ob231=-H32.'+X2-Mgse_tear.'*Bij_tear.';
                                    Ob241=-Dij_tear;
                                    Ob251=lambda*H22;
                                    Ob261=zeros(6,1);
                                    Ob271=-J_tear*X2;
                                    Ob281=kappa* Mgse_tear.'*Rbij.';
                                    
                                    Ob311=Ob131.';
                                    Ob321=Ob231.';
                                    Ob331=lambda*Q2+X2.'+X2;
                                    Ob341=-Dij_tear;
                                    Ob351=lambda*H32;
                                    Ob361=zeros(6,1);
                                    Ob371=-J_tear*X2;
                                    Ob381=zeros(6,6);
                                    
                                    Ob411=Ob141.';
                                    Ob421=Ob241.';
                                    Ob431=Ob341.';
                                    Ob441=-gamma*gamma*eye(1);
                                    Ob451=zeros(1,6);
                                    Ob461=zeros(1,1);
                                    Ob471=zeros(1,6);
                                    Ob481=zeros(1,6);
                                    
                                    Ob511=Ob151.';
                                    Ob521=Ob251.';
                                    Ob531=Ob351.';
                                    Ob541=Ob451.';
                                    Ob551=-lambda*Q2;
                                    Ob561=zeros(6,1);
                                    Ob571=zeros(6,6);
                                    Ob581=zeros(6,6);
                                    
                                    
                                    
                                    Ob611=Ob161.';
                                    Ob621=Ob261.';
                                    Ob631=Ob361.';
                                    Ob641=Ob461.';
                                    Ob651=Ob561.';
                                    Ob661=-eye(1);
                                    Ob671=zeros(1,6);
                                    Ob681=zeros(1,6);
                                    
                                    Ob711=Ob171.';
                                    Ob721=Ob271.';
                                    Ob731=Ob371.';
                                    Ob741=Ob471.';
                                    Ob751=Ob571.';
                                    Ob761=Ob671.';
                                    Ob771=-kappa*eye(6);
                                    Ob781=zeros(6,6);
                                    
                                    Ob811=Ob181.';
                                    Ob821=Ob281.';
                                    Ob831=Ob381.';
                                    Ob841=Ob481.';
                                    Ob851=Ob581.';
                                    Ob861=Ob681.';
                                    Ob871=Ob781.';
                                    Ob881=-kappa*eye(6);
                                    
                                    O21211=[Ob111  Ob121  Ob131  Ob141  Ob151  Ob161  Ob171  Ob181;
                                        Ob211  Ob221  Ob231  Ob241  Ob251  Ob261  Ob271  Ob281;
                                        Ob311  Ob321  Ob331  Ob341  Ob351  Ob361  Ob371  Ob381;
                                        Ob411  Ob421  Ob431  Ob441  Ob451  Ob461  Ob471  Ob481;
                                        Ob511  Ob521  Ob531  Ob541  Ob551  Ob561  Ob571  Ob581;
                                        Ob611  Ob621  Ob631  Ob641  Ob651  Ob661  Ob671  Ob681;
                                        Ob711  Ob721  Ob731  Ob741  Ob751  Ob761  Ob771  Ob781;
                                        Ob811  Ob821  Ob831  Ob841  Ob851  Ob861  Ob871  Ob881];
                                    
                                    
                                    %2121����
                                    SOS9=s3.'*(-[O12121]-[O21211]-eps3.*eye(38))*s3;
                                    % �ŧi�������x�}
                                    
                                    prog = sosineq(prog,SOS9,'sparsemultipartite',{ vars.' s3.'});
                                end
                                for subblock=8
                                    %% í�w����3-53,i<j=g=s,(4*3)x(4*3):
                                    %�X��2112�M1221
                                    % Aij_tear=inv(L1)*T*(mAi+nAj)*L1;
                                    % Bij_tear=inv(L1)*T*(mAi+nAj);
                                    % Gij_tear=mCi+nCj;
                                    % Dij_tear=inv(L1)*T*(mDi+nDj);
                                    % Mgse_tear=mMge+nMse;
                                    %2112�}�l
                                    Aij_tear=inv(L1)*T*(m*A2+n*A1)*L1;
                                    Bij_tear=inv(L1)*T*(m*B2+n*B1);
                                    Cij_tear=(m*C2+n*C1)*L1;
                                    Dij_tear=inv(L1)*T*(m*D2+n*D1);
                                    Mgse_tear=m*M12+n*M22;
                                    Sum=Sum22;
                                    
                                    Oa111=diff(P2,x2).*(A1(3,:)*x_h)+H12+H12.'-Aij_tear*X2-X2.'*Aij_tear+Sum;
                                    Oa121=-H12+H22.'-Bij_tear*Mgse_tear-X2.'*Aij_tear.';
                                    Oa131=P2+H32.'+X2-X2.'*Aij_tear;
                                    Oa141=-Dij_tear;
                                    Oa151=lambda*H12;
                                    Oa161=(Cij_tear*X2).';
                                    Oa171=-J_tear*X2;
                                    Oa181=kappa*X2.'*Raij.';
                                    
                                    Oa211=Oa121.';
                                    Oa221=-H22-H22.'-Bij_tear*Mgse_tear-Mgse_tear.'*Bij_tear.';
                                    Oa231=-H32.'+X2-Mgse_tear.'*Bij_tear.';
                                    Oa241=-Dij_tear;
                                    Oa251=lambda*H22;
                                    Oa261=zeros(6,1);
                                    Oa271=-J_tear*X2;
                                    Oa281=kappa* Mgse_tear.'*Rbij.';
                                    
                                    Oa311=Oa131.';
                                    Oa321=Oa231.';
                                    Oa331=lambda*Q2+X2.'+X2;
                                    Oa341=-Dij_tear;
                                    Oa351=lambda*H32;
                                    Oa361=zeros(6,1);
                                    Oa371=-J_tear*X2;
                                    Oa381=zeros(6,6);
                                    
                                    Oa411=Oa141.';
                                    Oa421=Oa241.';
                                    Oa431=Oa341.';
                                    Oa441=-gamma*gamma*eye(1);
                                    Oa451=zeros(1,6);
                                    Oa461=zeros(1,1);
                                    Oa471=zeros(1,6);
                                    Oa481=zeros(1,6);
                                    
                                    Oa511=Oa151.';
                                    Oa521=Oa251.';
                                    Oa531=Oa351.';
                                    Oa541=Oa451.';
                                    Oa551=-lambda*Q2;
                                    Oa561=zeros(6,1);
                                    Oa571=zeros(6,6);
                                    Oa581=zeros(6,6);
                                    
                                    
                                    
                                    Oa611=Oa161.';
                                    Oa621=Oa261.';
                                    Oa631=Oa361.';
                                    Oa641=Oa461.';
                                    Oa651=Oa561.';
                                    Oa661=-eye(1);
                                    Oa671=zeros(1,6);
                                    Oa681=zeros(1,6);
                                    
                                    Oa711=Oa171.';
                                    Oa721=Oa271.';
                                    Oa731=Oa371.';
                                    Oa741=Oa471.';
                                    Oa751=Oa571.';
                                    Oa761=Oa671.';
                                    Oa771=-kappa*eye(6);
                                    Oa781=zeros(6,6);
                                    
                                    Oa811=Oa181.';
                                    Oa821=Oa281.';
                                    Oa831=Oa381.';
                                    Oa841=Oa481.';
                                    Oa851=Oa581.';
                                    Oa861=Oa681.';
                                    Oa871=Oa781.';
                                    Oa881=-kappa*eye(6);
                                    
                                    O21121=[Oa111  Oa121  Oa131  Oa141  Oa151  Oa161  Oa171  Oa181;
                                        Oa211  Oa221  Oa231  Oa241  Oa251  Oa261  Oa271  Oa281;
                                        Oa311  Oa321  Oa331  Oa341  Oa351  Oa361  Oa371  Oa381;
                                        Oa411  Oa421  Oa431  Oa441  Oa451  Oa461  Oa471  Oa481;
                                        Oa511  Oa521  Oa531  Oa541  Oa551  Oa561  Oa571  Oa581;
                                        Oa611  Oa621  Oa631  Oa641  Oa651  Oa661  Oa671  Oa681;
                                        Oa711  Oa721  Oa731  Oa741  Oa751  Oa761  Oa771  Oa781;
                                        Oa811  Oa821  Oa831  Oa841  Oa851  Oa861  Oa871  Oa881];
                                    
                                    %2112����
                                    %1221�}�l
                                    Aij_tear=inv(L1)*T*(m*A1+n*A2)*L1;
                                    Bij_tear=inv(L1)*T*(m*B1+n*B2);
                                    Cij_tear=(m*C1+n*C2)*L1;
                                    Dij_tear=inv(L1)*T*(m*D1+n*D2);
                                    Mgse_tear=m*M22+n*M12;
                                    Sum=Sum12;
                                    
                                    Ob111=diff(P2,x2).*(A1(3,:)*x_h)+H12+H12.'-Aij_tear*X2-X2.'*Aij_tear+Sum;
                                    Ob121=-H12+H22.'-Bij_tear*Mgse_tear-X2.'*Aij_tear.';
                                    Ob131=P2+H32.'+X2-X2.'*Aij_tear;
                                    Ob141=-Dij_tear;
                                    Ob151=lambda*H12;
                                    Ob161=(Cij_tear*X2).';
                                    Ob171=-J_tear*X2;
                                    Ob181=kappa*X2.'*Raij.';
                                    
                                    Ob211=Ob121.';
                                    Ob221=-H22-H22.'-Bij_tear*Mgse_tear-Mgse_tear.'*Bij_tear.';
                                    Ob231=-H32.'+X2-Mgse_tear.'*Bij_tear.';
                                    Ob241=-Dij_tear;
                                    Ob251=lambda*H22;
                                    Ob261=zeros(6,1);
                                    Ob271=-J_tear*X2;
                                    Ob281=kappa* Mgse_tear.'*Rbij.';
                                    
                                    Ob311=Ob131.';
                                    Ob321=Ob231.';
                                    Ob331=lambda*Q2+X2.'+X2;
                                    Ob341=-Dij_tear;
                                    Ob351=lambda*H32;
                                    Ob361=zeros(6,1);
                                    Ob371=-J_tear*X2;
                                    Ob381=zeros(6,6);
                                    
                                    Ob411=Ob141.';
                                    Ob421=Ob241.';
                                    Ob431=Ob341.';
                                    Ob441=-gamma*gamma*eye(1);
                                    Ob451=zeros(1,6);
                                    Ob461=zeros(1,1);
                                    Ob471=zeros(1,6);
                                    Ob481=zeros(1,6);
                                    
                                    Ob511=Ob151.';
                                    Ob521=Ob251.';
                                    Ob531=Ob351.';
                                    Ob541=Ob451.';
                                    Ob551=-lambda*Q2;
                                    Ob561=zeros(6,1);
                                    Ob571=zeros(6,6);
                                    Ob581=zeros(6,6);
                                    
                                    
                                    
                                    Ob611=Ob161.';
                                    Ob621=Ob261.';
                                    Ob631=Ob361.';
                                    Ob641=Ob461.';
                                    Ob651=Ob561.';
                                    Ob661=-eye(1);
                                    Ob671=zeros(1,6);
                                    Ob681=zeros(1,6);
                                    
                                    Ob711=Ob171.';
                                    Ob721=Ob271.';
                                    Ob731=Ob371.';
                                    Ob741=Ob471.';
                                    Ob751=Ob571.';
                                    Ob761=Ob671.';
                                    Ob771=-kappa*eye(6);
                                    Ob781=zeros(6,6);
                                    
                                    Ob811=Ob181.';
                                    Ob821=Ob281.';
                                    Ob831=Ob381.';
                                    Ob841=Ob481.';
                                    Ob851=Ob581.';
                                    Ob861=Ob681.';
                                    Ob871=Ob781.';
                                    Ob881=-kappa*eye(6);
                                    
                                    O12211=[Ob111  Ob121  Ob131  Ob141  Ob151  Ob161  Ob171  Ob181;
                                        Ob211  Ob221  Ob231  Ob241  Ob251  Ob261  Ob271  Ob281;
                                        Ob311  Ob321  Ob331  Ob341  Ob351  Ob361  Ob371  Ob381;
                                        Ob411  Ob421  Ob431  Ob441  Ob451  Ob461  Ob471  Ob481;
                                        Ob511  Ob521  Ob531  Ob541  Ob551  Ob561  Ob571  Ob581;
                                        Ob611  Ob621  Ob631  Ob641  Ob651  Ob661  Ob671  Ob681;
                                        Ob711  Ob721  Ob731  Ob741  Ob751  Ob761  Ob771  Ob781;
                                        Ob811  Ob821  Ob831  Ob841  Ob851  Ob861  Ob871  Ob881];
                                    
                                    
                                    %1221����
                                    SOS10=s3.'*(-[O21121]-[O12211]-eps3.*eye(38))*s3;
                                    % �ŧi�������x�}
                                    
                                    prog = sosineq(prog,SOS10,'sparsemultipartite',{ vars.' s3.'});
                                end
                                for subblock=9
                                    %% í�w����3-53,i<j=g=s,(4*3)x(4*3):
                                    %�X��2212�M1121
                                    % Aij_tear=inv(L1)*T*(mAi+nAj)*L1;
                                    % Bij_tear=inv(L1)*T*(mAi+nAj);
                                    % Gij_tear=mCi+nCj;
                                    % Dij_tear=inv(L1)*T*(mDi+nDj);
                                    % Mgse_tear=mMge+nMse;
                                    %2212�}�l
                                    Aij_tear=inv(L1)*T*(m*A2+n*A2)*L1;
                                    Bij_tear=inv(L1)*T*(m*B2+n*B2);
                                    Cij_tear=(m*C2+n*C2)*L1;
                                    Dij_tear=inv(L1)*T*(m*D2+n*D2);
                                    Mgse_tear=m*M12+n*M22;
                                    Sum=Sum22;
                                    
                                    Oa111=diff(P2,x2).*(A1(3,:)*x_h)+H12+H12.'-Aij_tear*X2-X2.'*Aij_tear+Sum;
                                    Oa121=-H12+H22.'-Bij_tear*Mgse_tear-X2.'*Aij_tear.';
                                    Oa131=P2+H32.'+X2-X2.'*Aij_tear;
                                    Oa141=-Dij_tear;
                                    Oa151=lambda*H12;
                                    Oa161=(Cij_tear*X2).';
                                    Oa171=-J_tear*X2;
                                    Oa181=kappa*X2.'*Raij.';
                                    
                                    Oa211=Oa121.';
                                    Oa221=-H22-H22.'-Bij_tear*Mgse_tear-Mgse_tear.'*Bij_tear.';
                                    Oa231=-H32.'+X2-Mgse_tear.'*Bij_tear.';
                                    Oa241=-Dij_tear;
                                    Oa251=lambda*H22;
                                    Oa261=zeros(6,1);
                                    Oa271=-J_tear*X2;
                                    Oa281=kappa* Mgse_tear.'*Rbij.';
                                    
                                    Oa311=Oa131.';
                                    Oa321=Oa231.';
                                    Oa331=lambda*Q2+X2.'+X2;
                                    Oa341=-Dij_tear;
                                    Oa351=lambda*H32;
                                    Oa361=zeros(6,1);
                                    Oa371=-J_tear*X2;
                                    Oa381=zeros(6,6);
                                    
                                    Oa411=Oa141.';
                                    Oa421=Oa241.';
                                    Oa431=Oa341.';
                                    Oa441=-gamma*gamma*eye(1);
                                    Oa451=zeros(1,6);
                                    Oa461=zeros(1,1);
                                    Oa471=zeros(1,6);
                                    Oa481=zeros(1,6);
                                    
                                    Oa511=Oa151.';
                                    Oa521=Oa251.';
                                    Oa531=Oa351.';
                                    Oa541=Oa451.';
                                    Oa551=-lambda*Q2;
                                    Oa561=zeros(6,1);
                                    Oa571=zeros(6,6);
                                    Oa581=zeros(6,6);
                                    
                                    
                                    
                                    Oa611=Oa161.';
                                    Oa621=Oa261.';
                                    Oa631=Oa361.';
                                    Oa641=Oa461.';
                                    Oa651=Oa561.';
                                    Oa661=-eye(1);
                                    Oa671=zeros(1,6);
                                    Oa681=zeros(1,6);
                                    
                                    Oa711=Oa171.';
                                    Oa721=Oa271.';
                                    Oa731=Oa371.';
                                    Oa741=Oa471.';
                                    Oa751=Oa571.';
                                    Oa761=Oa671.';
                                    Oa771=-kappa*eye(6);
                                    Oa781=zeros(6,6);
                                    
                                    Oa811=Oa181.';
                                    Oa821=Oa281.';
                                    Oa831=Oa381.';
                                    Oa841=Oa481.';
                                    Oa851=Oa581.';
                                    Oa861=Oa681.';
                                    Oa871=Oa781.';
                                    Oa881=-kappa*eye(6);
                                    
                                    O22121=[Oa111  Oa121  Oa131  Oa141  Oa151  Oa161  Oa171  Oa181;
                                        Oa211  Oa221  Oa231  Oa241  Oa251  Oa261  Oa271  Oa281;
                                        Oa311  Oa321  Oa331  Oa341  Oa351  Oa361  Oa371  Oa381;
                                        Oa411  Oa421  Oa431  Oa441  Oa451  Oa461  Oa471  Oa481;
                                        Oa511  Oa521  Oa531  Oa541  Oa551  Oa561  Oa571  Oa581;
                                        Oa611  Oa621  Oa631  Oa641  Oa651  Oa661  Oa671  Oa681;
                                        Oa711  Oa721  Oa731  Oa741  Oa751  Oa761  Oa771  Oa781;
                                        Oa811  Oa821  Oa831  Oa841  Oa851  Oa861  Oa871  Oa881];
                                    
                                    %2212����
                                    %1121�}�l
                                    Aij_tear=inv(L1)*T*(m*A1+n*A1)*L1;
                                    Bij_tear=inv(L1)*T*(m*B1+n*B1);
                                    Cij_tear=(m*C1+n*C1)*L1;
                                    Dij_tear=inv(L1)*T*(m*D1+n*D1);
                                    Mgse_tear=m*M22+n*M12;
                                    Sum=Sum12;
                                    
                                    Ob111=diff(P2,x2).*(A1(3,:)*x_h)+H12+H12.'-Aij_tear*X2-X2.'*Aij_tear+Sum;
                                    Ob121=-H12+H22.'-Bij_tear*Mgse_tear-X2.'*Aij_tear.';
                                    Ob131=P2+H32.'+X2-X2.'*Aij_tear;
                                    Ob141=-Dij_tear;
                                    Ob151=lambda*H12;
                                    Ob161=(Cij_tear*X2).';
                                    Ob171=-J_tear*X2;
                                    Ob181=kappa*X2.'*Raij.';
                                    
                                    Ob211=Ob121.';
                                    Ob221=-H22-H22.'-Bij_tear*Mgse_tear-Mgse_tear.'*Bij_tear.';
                                    Ob231=-H32.'+X2-Mgse_tear.'*Bij_tear.';
                                    Ob241=-Dij_tear;
                                    Ob251=lambda*H22;
                                    Ob261=zeros(6,1);
                                    Ob271=-J_tear*X2;
                                    Ob281=kappa* Mgse_tear.'*Rbij.';
                                    
                                    Ob311=Ob131.';
                                    Ob321=Ob231.';
                                    Ob331=lambda*Q2+X2.'+X2;
                                    Ob341=-Dij_tear;
                                    Ob351=lambda*H32;
                                    Ob361=zeros(6,1);
                                    Ob371=-J_tear*X2;
                                    Ob381=zeros(6,6);
                                    
                                    Ob411=Ob141.';
                                    Ob421=Ob241.';
                                    Ob431=Ob341.';
                                    Ob441=-gamma*gamma*eye(1);
                                    Ob451=zeros(1,6);
                                    Ob461=zeros(1,1);
                                    Ob471=zeros(1,6);
                                    Ob481=zeros(1,6);
                                    
                                    Ob511=Ob151.';
                                    Ob521=Ob251.';
                                    Ob531=Ob351.';
                                    Ob541=Ob451.';
                                    Ob551=-lambda*Q2;
                                    Ob561=zeros(6,1);
                                    Ob571=zeros(6,6);
                                    Ob581=zeros(6,6);
                                    
                                    
                                    
                                    Ob611=Ob161.';
                                    Ob621=Ob261.';
                                    Ob631=Ob361.';
                                    Ob641=Ob461.';
                                    Ob651=Ob561.';
                                    Ob661=-eye(1);
                                    Ob671=zeros(1,6);
                                    Ob681=zeros(1,6);
                                    
                                    Ob711=Ob171.';
                                    Ob721=Ob271.';
                                    Ob731=Ob371.';
                                    Ob741=Ob471.';
                                    Ob751=Ob571.';
                                    Ob761=Ob671.';
                                    Ob771=-kappa*eye(6);
                                    Ob781=zeros(6,6);
                                    
                                    Ob811=Ob181.';
                                    Ob821=Ob281.';
                                    Ob831=Ob381.';
                                    Ob841=Ob481.';
                                    Ob851=Ob581.';
                                    Ob861=Ob681.';
                                    Ob871=Ob781.';
                                    Ob881=-kappa*eye(6);
                                    
                                    O11211=[Ob111  Ob121  Ob131  Ob141  Ob151  Ob161  Ob171  Ob181;
                                        Ob211  Ob221  Ob231  Ob241  Ob251  Ob261  Ob271  Ob281;
                                        Ob311  Ob321  Ob331  Ob341  Ob351  Ob361  Ob371  Ob381;
                                        Ob411  Ob421  Ob431  Ob441  Ob451  Ob461  Ob471  Ob481;
                                        Ob511  Ob521  Ob531  Ob541  Ob551  Ob561  Ob571  Ob581;
                                        Ob611  Ob621  Ob631  Ob641  Ob651  Ob661  Ob671  Ob681;
                                        Ob711  Ob721  Ob731  Ob741  Ob751  Ob761  Ob771  Ob781;
                                        Ob811  Ob821  Ob831  Ob841  Ob851  Ob861  Ob871  Ob881];
                                    
                                    %2211����
                                    SOS11=s3.'*(-[O22121]-[O11211]-eps3.*eye(38))*s3;
                                    % �ŧi�������x�}
                                    
                                    prog = sosineq(prog,SOS11,'sparsemultipartite',{ vars.' s3.'});
                                end
                            end
                        end
                        
                        
                        %}
                        %% �D�ѯx�}
                        for block=5
                            % �I�sSOS�D�Ѿ�
                            %N=1
                            prog=sossolve(prog);
                            % �o�춷�D�Ѫ��ܼ�
                            P1=sosgetsol(prog,P1)
                            X1=sosgetsol(prog,X1)
                            Q1=sosgetsol(prog,Q1)
                            H11=sosgetsol(prog,H11)
                            H21=sosgetsol(prog,H21)
                            H31=sosgetsol(prog,H31)
                            M11=sosgetsol(prog,M11)
                            M21=sosgetsol(prog,M21)
                            F11=M11*inv(X1)
                            F21=M21*inv(X1)
                            eig_P1=eig(P1)
                            eig_Q1=eig(Q1)
                            %N=2
                            % �o�춷�D�Ѫ��ܼ�
                            P2=sosgetsol(prog,P2)
                            X2=sosgetsol(prog,X2)
                            Q2=sosgetsol(prog,Q2)
                            H12=sosgetsol(prog,H12)
                            H22=sosgetsol(prog,H22)
                            H32=sosgetsol(prog,H32)
                            M12=sosgetsol(prog,M12)
                            M22=sosgetsol(prog,M22)
                            F12=M12*inv(X2)
                            F22=M22*inv(X2)
                            eig_P2=eig(P2)
                            eig_Q2=eig(Q2)
                        end
                        %% �x�s�S�x�Ȥj��0����
                        for block=6
                            if N == 1
                                if M11(1,1)~=0 || M11(1,2)~=0 || M11(1,3)~=0 ||M11(1,4)~=0 || M11(1,5)~=0 || M11(1,6)~=0 || M11(2,1)~=0 || M11(2,2)~=0 || M11(2,3)~=0 || M11(2,4)~=0 || M11(2,5)~=0 || M11(2,6)~=0
                                    if M21(1,1)~=0 || M21(1,2)~=0 || M21(1,3)~=0 ||M21(1,4)~=0 || M21(1,5)~=0 || M21(1,6)~=0 || M21(2,1)~=0 || M21(2,2)~=0 || M21(2,3)~=0 || M21(2,4)~=0 || M21(2,5)~=0 || M21(2,6)~=0
                                        
                                        %�ΨӰO��K���x�}
                                        count_1 = count_1+1;
                                        %�S�x�ȥ����~�O����
                                        if eig_P1(1,1)>0 && eig_P1(2,1)>0 && eig_P1(3,1)>0 && eig_Q1(1,1)>0 && eig_Q1(2,1)>0 && eig_Q1(3,1)>0
                                            
                                            sol_space(count_1,1)=m;
                                            sol_space(count_1,2)=n;
                                            sol_space(count_1,3)=lambda;
                                            sol_space(count_1,4)=gamma;
                                            sol_space(count_1,5)=kappa;
                                            
                                            
                                            %K11,K21�s�b�x�}��
                                            fin_K1_data(count_1*2-1,1)=F11(1,1);   fin_K1_data(count_1*2-1,4)=F21(1,1);
                                            fin_K1_data(count_1*2-1,2)=F11(1,2);   fin_K1_data(count_1*2-1,5)=F21(1,2);
                                            fin_K1_data(count_1*2-1,3)=F11(1,3);   fin_K1_data(count_1*2-1,6)=F21(1,3);
                                            fin_K1_data(count_1*2,1)=F11(2,1); fin_K1_data(count_1*2,4)=F21(2,1);
                                            fin_K1_data(count_1*2,2)=F11(2,2); fin_K1_data(count_1*2,5)=F21(2,2);
                                            fin_K1_data(count_1*2,3)=F11(2,3); fin_K1_data(count_1*2,6)=F21(2,3);
                                        end
                                    end
                                end
                            end
                            if N == 2
                                if M11(1,1)~=0 || M11(1,2)~=0 || M11(1,3)~=0 ||M11(1,4)~=0 || M11(1,5)~=0 || M11(1,6)~=0 || M11(2,1)~=0 || M11(2,2)~=0 || M11(2,3)~=0 || M11(2,4)~=0 || M11(2,5)~=0 || M11(2,6)~=0
                                    if M21(1,1)~=0 || M21(1,2)~=0 || M21(1,3)~=0 ||M21(1,4)~=0 || M21(1,5)~=0 || M21(1,6)~=0 || M21(2,1)~=0 || M21(2,2)~=0 || M21(2,3)~=0 || M21(2,4)~=0 || M21(2,5)~=0 || M21(2,6)~=0
                                        if M12(1,1)~=0 || M12(1,2)~=0 || M12(1,3)~=0 ||M12(1,4)~=0 || M12(1,5)~=0 || M12(1,6)~=0 || M12(2,1)~=0 || M12(2,2)~=0 || M12(2,3)~=0 || M12(2,4)~=0 || M12(2,5)~=0 || M12(2,6)~=0
                                            if M22(1,1)~=0 || M22(1,2)~=0 || M22(1,3)~=0 ||M22(1,4)~=0 || M22(1,5)~=0 || M22(1,6)~=0 || M22(2,1)~=0 || M22(2,2)~=0 || M22(2,3)~=0 || M22(2,4)~=0 || M22(2,5)~=0 || M22(2,6)~=0
                                                %�ΨӰO��K���x�}
                                                count_1 = count_1+1;
                                                %�S�x�ȥ����~�O����
                                                if eig_P1(1,1)>0 && eig_P1(2,1)>0 && eig_P1(3,1)>0 && eig_Q1(1,1)>0 && eig_Q1(2,1)>0 && eig_Q1(3,1)>0
                                                    
                                                    sol_space(count_1,1)=m;
                                                    sol_space(count_1,2)=n;
                                                    sol_space(count_1,3)=lambda;
                                                    sol_space(count_1,4)=gamma;
                                                    sol_space(count_1,5)=kappa;
                                                    
                                                    
                                                    %K11,K21�s�b�x�}��
                                                    fin_K1_data(count_1*2-1,1)=F11(1,1);   fin_K1_data(count_1*2-1,1)=F21(1,1);
                                                    fin_K1_data(count_1*2-1,2)=F11(1,2);   fin_K1_data(count_1*2-1,2)=F21(1,2);
                                                    fin_K1_data(count_1*2-1,3)=F11(1,3);   fin_K1_data(count_1*2-1,3)=F21(1,3);
                                                    fin_K1_data(count_1*2-1,4)=F11(1,4);   fin_K1_data(count_1*2-1,4)=F21(1,4);
                                                    fin_K1_data(count_1*2-1,5)=F11(1,5);   fin_K1_data(count_1*2-1,5)=F21(1,5);
                                                    fin_K1_data(count_1*2-1,6)=F11(1,6);   fin_K1_data(count_1*2-1,6)=F21(1,6);
                                                    fin_K1_data(count_1*2,1)=F11(2,1); fin_K1_data(count_1*2,1)=F21(2,1);
                                                    fin_K1_data(count_1*2,2)=F11(2,2); fin_K1_data(count_1*2,2)=F21(2,2);
                                                    fin_K1_data(count_1*2,3)=F11(2,3); fin_K1_data(count_1*2,3)=F21(2,3);
                                                    fin_K1_data(count_1*2,4)=F11(2,4); fin_K1_data(count_1*2,4)=F21(2,4);
                                                    fin_K1_data(count_1*2,5)=F11(2,5); fin_K1_data(count_1*2,5)=F21(2,5);
                                                    fin_K1_data(count_1*2,6)=F11(2,6); fin_K1_data(count_1*2,6)=F21(2,6);
                                                    %K12,K22�s�b�x�}��
                                                    fin_K2_data(count_1*2-1,1)=F12(1,1);   fin_K2_data(count_1*2-1,4)=F22(1,1);
                                                    fin_K2_data(count_1*2-1,2)=F12(1,2);   fin_K2_data(count_1*2-1,5)=F22(1,2);
                                                    fin_K2_data(count_1*2-1,3)=F12(1,3);   fin_K2_data(count_1*2-1,6)=F22(1,3);
                                                    fin_K2_data(count_1*2-1,4)=F12(1,1);   fin_K2_data(count_1*2-1,4)=F22(1,1);
                                                    fin_K2_data(count_1*2-1,5)=F12(1,2);   fin_K2_data(count_1*2-1,5)=F22(1,2);
                                                    fin_K2_data(count_1*2-1,6)=F12(1,3);   fin_K2_data(count_1*2-1,6)=F22(1,3);
                                                    fin_K2_data(count_1*2,1)=F12(2,1); fin_K2_data(count_1*2,4)=F22(2,1);
                                                    fin_K2_data(count_1*2,2)=F12(2,2); fin_K2_data(count_1*2,5)=F22(2,2);
                                                    fin_K2_data(count_1*2,3)=F12(2,3); fin_K2_data(count_1*2,6)=F22(2,3);
                                                    fin_K2_data(count_1*2,4)=F12(2,1); fin_K2_data(count_1*2,4)=F22(2,1);
                                                    fin_K2_data(count_1*2,5)=F12(2,2); fin_K2_data(count_1*2,5)=F22(2,2);
                                                    fin_K2_data(count_1*2,6)=F12(2,3); fin_K2_data(count_1*2,6)=F22(2,3);
                                                    
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end




toc