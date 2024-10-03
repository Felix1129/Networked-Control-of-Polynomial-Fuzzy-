# Networked Control of Polynomial Fuzzy Systems with Time-Delay Based on Piecewise Lyapunov Functions
## Summary
The networked control systems (NCSs) of polynomial fuzzy time-delay systems are based on piecewise polynomial Lyapunov functions (PPLF).
Applied to the path-tracking control of mobile robots and quadrotors. 

The number of fuzzy rules could be reduced to two by using the novel modeling method, which greatly reduces the hardware cost. 
Both mobile robots and quadrotors had the problem of time delay in the plant. 
The time delay system developed in this work could tolerate the time delay of the plant through the modeling method of calculating the time delay of the plant additionally. 
The Lyapunov-Krasovskii piecewise polynomial functional function was designed using the minimum-type polynomial Lyapunov piecewise function. 
By reducing the stable condition number, the gain solution space could be expanded.

Based on the piecewise polynomial Lyapunov-Krasovskii function, considering the time delay and package dropout caused by the network, a robust stability condition is proposed in this work with calculated time delay in the SOS-based relaxation stabilization condition, and an SOS-based formation relaxation condition. 
The theorems proposed in this work consider the network-caused package dropout and time delay. 


## *Networked Control Systems*
### The block diagram of NCS
![image](https://github.com/user-attachments/assets/5add960e-4817-4cfd-8177-cb49e7b43dca)

Let `tk, tk+1, …, tk+q (k=0,1,2,…) `indicate the sampling instants, 
`hp=tk+1 − tk` means the sampling period `(q=1,2,…,qmax)` represent the total time delay caused by the network.
Considering the delay and package dropout caused by the 
network, two parameters are often used to design the NCS 
controller. The first is the maximum allowable delay bound `ρ`, 
which is the maximum permissible control interval from the 
moment the sensors transmit the sampled data measured by the 
plant to the moment the actuators send the control signal. The 
second parameter is the maximum allowable transfer interval `δ`, 
which is a period of time; in other words, if a control signal is 
transmitted at a point in time` k k<sub>t</sub>`. Based on the 
two definitions above,` λ = ρ + δ` can be defined as the maximum 
allowable control interval to address both network-induced 
delay and package dropout problems simultaneously. 

### Control signal time-sequence diagram
![image](https://github.com/user-attachments/assets/58db72ed-4840-4c50-9f31-8b72a033f769)

For a given `λ`, inequality is assumed as follows: 

![image](https://github.com/user-attachments/assets/d590c72a-651e-4a9d-a7ff-4d2976b9f3ee)

where `λ` is an upper bound that ensures the closed-loop system maintains stability. 
Moreover, the following is the definition of the number of data: 

![image](https://github.com/user-attachments/assets/c3e0dd5a-b805-499c-aef8-0deb0033980e)


## *Minimum-Type Piecewise Polynomial Lyapunov Function*
 The minimum-type PPLF is obtained as follows: 
 
 ![image](https://github.com/user-attachments/assets/928b3d59-4404-4f92-aca3-40227e55c17f)
 

## *Novel Output Feedback Piecewise Polynomial Fuzzy Networked Control Time-Delay Systems*
closed-loop novel piecewise polynomial fuzzy networked control time-delay system, 
combining a novel polynomial fuzzy model with external 
disturbances and model uncertainties and a novel piecewise controller based on output feedback.

### *Model Rule i*
![image](https://github.com/user-attachments/assets/3b8642fa-d2cc-463c-8427-c77323d4df1d)

where `J , R<sub>ai</sub> , and R<sub>bi</sub>` are constant matrices with 
corresponding dimension. is a time-varying matrices bounded 
by `K<sup>T</sup><sub>(t)</sub>K<sub>(t)</sub><I`. 
The defuzzification result of novel polynomial fuzzy NCSs can be indicated as
![image](https://github.com/user-attachments/assets/eae33f61-dbf5-4074-b7cd-2393104ab09a)

### *Control Rule i*
The defuzzification of the control rules is shown as follows:
![image](https://github.com/user-attachments/assets/5601d710-cc31-4450-a048-e86184fb9c4b)

The transformed output feedback system of the defuzzification result of novel polynomial fuzzy NCSs can be obtained as
![image](https://github.com/user-attachments/assets/bf9e815e-bdcf-4272-8f13-e52d5b3d9494)

# [Simulation](https://github.com/Felix1129/Networked-Control-of-Polynomial-Fuzzy-/tree/main/Simulation)

## [Package dropout](https://github.com/Felix1129/Networked-Control-of-Polynomial-Fuzzy-/tree/main/Simulation/Packet%20dropout)
To prove that the negative impact of package dropout on 
NCSs and package dropout is stochastic, we set the package 
`dropout rate  = 50%` and the disturbance to occur at 20 to 21 
s, with a pulse with an amplitude of 0.8. Since package 
dropout is a random event, the error time response graphs of `x`, 
`y`, and `θ` are presented in four simulations. The red line is a 
discrete signal and the blue line is a continuous signal. Based 
on these, an apparent simulation diagram is provided that the 
proposed new piecewise polynomial fuzzy controller can still 
correctly track continuous signals when the network causes 
package dropout. This verifies that the proposed novel 
piecewise polynomial fuzzy controller can effectively control 
the network system while ensuring the stability of the system.
![image](https://github.com/user-attachments/assets/9b91badd-3847-4154-ac3a-44ab70082afa)
![image](https://github.com/user-attachments/assets/d13e4ba0-2918-418c-9966-b2edcd791855)
![image](https://github.com/user-attachments/assets/617ca08e-7351-4011-87d1-3318b5637ddf)

## [Feasible solution space](https://github.com/Felix1129/Networked-Control-of-Polynomial-Fuzzy-/tree/main/Simulation/Packet%20dropout)
A comparison of the proposed method with the existing 
method indicated an increase in both the number of fuzzy rules 
and stability conditions. It can be verified by computing time in 
Table 1. 

![image](https://github.com/user-attachments/assets/2b72db25-8dfd-453b-a40d-9d70333a4f71)

![image](https://github.com/user-attachments/assets/66396be3-4ea3-42fe-a3b8-9be1972b6556)

## [*Mobile Robot*](https://github.com/Felix1129/Networked-Control-of-Polynomial-Fuzzy-/tree/main/Simulation/Mobile%20Robot)
The error of the proposed method is shown in red, 
while the error of the existing method is shown in blue. The segment of the Lyapunov 
function is equal to 2. In the figure, the red line shows a faster 
response speed at the beginning than the blue line does, and the 
magnitude of oscillation is smaller. The tracking 
trajectory of the mobile robot in the X-Y plane. The black line 
represents the desired trajectory, the red line means the 
proposed method tracking trajectory, and the blue line means 
the existing method tracking trajectory. According to the figure, 
the red line demonstrates better tracking performance than the 
blue line. 
![image](https://github.com/user-attachments/assets/ed9001e5-2121-45bb-b737-65e669623d5f)
![image](https://github.com/user-attachments/assets/1e7d5f2f-1143-4989-aeca-77390f2f3cb7)

## [*Quadrotor*](https://github.com/Felix1129/Networked-Control-of-Polynomial-Fuzzy-/tree/main/Simulation/Quadrotor)
The red line is the proposed novel 
piecewise polynomial fuzzy time-delay network controller, and 
the blue line is the existing piecewise polynomial fuzzy time delay network controller. 
In the figure, the response and convergence speed of the novel piecewise polynomial fuzzy 
delay controller is better than that of the existing delay controller, 
and it also converges faster than the existing controller when external disturbances. 
The integrated square error of the coordinate and attitude of the 
proposed methods are smaller than those of the existing ones.
![image](https://github.com/user-attachments/assets/8498e59d-f05f-43e2-ac9f-d6e5a4e276ea)
![image](https://github.com/user-attachments/assets/bb485237-c7c5-4971-9f8b-7a306a5801ed)
![image](https://github.com/user-attachments/assets/d63f070c-bb82-43b0-b3d4-843dbde5c35b)
![image](https://github.com/user-attachments/assets/7e6c0fd4-8330-4940-aded-e45af3347f88)



# [Experiment Results](https://github.com/Felix1129/Networked-Control-of-Polynomial-Fuzzy-/tree/main/Experiment%20Results)

## [*Mobile Robot*](https://github.com/Felix1129/Networked-Control-of-Polynomial-Fuzzy-/tree/main/Experiment%20Results/Mobile%20robot)
![image](https://github.com/user-attachments/assets/c79a2d6d-2c43-47e6-9ec4-ad0d2d9768c4)
![image](https://github.com/user-attachments/assets/fde7db6d-304d-4584-bd69-352c4d4f58f2)

## [*Quadrotor*](https://github.com/Felix1129/Networked-Control-of-Polynomial-Fuzzy-/tree/main/Experiment%20Results/Quadrotor)
![image](https://github.com/user-attachments/assets/e891b3c3-5b0c-40d4-a02a-84168150b18d)
![image](https://github.com/user-attachments/assets/ecff29fc-6972-491d-9999-92f00de0843d)
![image](https://github.com/user-attachments/assets/de0f9437-21eb-4102-bfc6-0162ecbe18a3)



