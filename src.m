%%%Final Project ENGR 313
%%%Spring '17
%%%Amr Adel Gouhar   900153482
%%%Eslam Ehab        900153480

clear %Clear all previous values in MATLAB
clc  % Clearing the screen


nodes = input('Please insert the number of nodes in your system '); %inserting the number of nodes in the network
streets = input('Please insert the number of streets in your system '); %inserting the number of streets in the network
Constraints = input('Please insert the maximum number of cars per street in order [x1 x2 x3 .. etc] in the following form:[2 3 4 5 ... etc]:\n'); %Taking constraints from the user %% if there is no constraint insert a very huge number
R = zeros(nodes,streets+1); %initializing a Matrix with zeros with size nodes x 

%Looping over the nodes to take the data from the user
for u = 1:nodes
   if u ==1 % Give instructions of how to take the unknowns from the user in the first time only
   fprintf('Please insert the indices of the unknowns in the following manner: x1(into the node) , x5(out of the node) , x6(into the node) is represented in the form: [1 -5 6]\n');
   end
   indices = input('Please insert the indices: '); %Take the indices from the user 
   if u == 1 % Give instructions of how to take the given values from the user in the first time only
   fprintf('Please insert the number of cars in each node in the following manner: 500(into the node) , 200(out of the node) , 300(into the node) is represented in the form: [500 -200 300]\n'); 
   end
   values = input('Please insert the number of cars: '); %Take the values from the user
%Inserting the Data into the intialized Matrix of zeros R;
   for t = 1:length(indices)
       num = abs(indices(t));
       if indices(t) > 0
          R(u,num) = 1;
       else R(u,num) = -1;
       end
   end
   R(u,streets+1) = - sum(values);
end
%Putting the data of R into x to solve it maintaing the Original Data in R
%to be used later if needed
x = R;
%n is the number of rows, m is the number of columns 
[n,m] = size(x);

%Looping over the rows to start solving the equations
for j=1:n-1
%Looping over the equations to pivot
    for r=2:n
        if x(j,j)==0
            tmp=x(1,:);
            x(1,:)=x(r,:);
            x(r,:)=tmp;
        end
    end
%Looping over the equations to solve the Down Traingle
    for i=j+1:n
        x(i,:)=x(i,:)-x(j,:)*(x(i,j)/x(j,j));
    end
end
%Looping back the rows to solve the UP Traingle
for j=n-1:-1:2
    for i=j-1:-1:1
        x(i,:)=x(i,:)-x(j,:)*(x(i,j)/x(j,j));
    end
end
%Preparing the Coefecients
for k=1:n-1
    x(k,:)=x(k,:)/x(k,k);
    P(k)=x(k,m);
end


flag = ones(1, streets); %Setting a vector of size streets into 1s

%If flag(3) = 1; this means that x3 weren't solved yet within the
%constraints, else: the program found an answer for x3 within the given
%constraints

zflag = sum(flag);  %Summing the vector of flags, if zflag = 0 then an answer is found for all the Unknowns within the constraints

counto = 0;  %this is is a counter to count a million iteration before deciding that the system has only one solution within the given constraints or doesn't have any

% RE1 is an array with the first solution
% setting an initail assumption for RE1 to 1
 for q = streets:-1:(nodes)
    RE1(q) = 1; 
 end
 
 % A loop to find the first set of answers to the network
while ((zflag ~= 0) && (counto < 10000000))  % it loops over until an answer is found OR the counto has counted a Million iteration and no answer is found

    % Looping over the unknowns to be assumed and find a random assumption
    for q = streets:-1:(nodes)
        tmpo = RE1(q);  % a temporary storage for the assumed result
        RE1(q) = randi(Constraints(q)); % a new assumption is made in the range between 1 and the given constraint
        
        % A loop is started to make sure that the new a assumption is
        % greater than the old assumption
        
        while (tmpo >= RE1(q) && counto <5000000)    % while the new assumption is less than the old assumption and the iteration counter is less than 500 Thousand
            TRAIL = randi(Constraints(q), [50,1]);  % Assume an array of valid assumptions with size 50
            if(max(TRAIL)+1 <= Constraints(q))  %% If the (maximum value + 1)is less than or equal the constraint 
                RE1(q) = max(TRAIL)+1;  % put the new maximum value +1 as an assumed answer
            else
                RE1(q) = max(TRAIL);    % put the new maximum value only as an assumed answer
            end

            flag(q) = 0; %set a flag for the assumed value to 0 %%(an answer is found for this unknown)
            counto = counto + 1; %increment the iteration counter
        end

        if counto >= 5000000 %if the iteration counter has exceeded 500000 then decrease the assumed value 
                        TRAIL = randi(RE1(q), [50,1]);  % assume an array of size 50 with valid assumptions that are less than the previous assumption
                        if(max(TRAIL)+1 <= Constraints(q)) %% If the (maximum value + 1)is less than or equal the constraint 
                            RE1(q) = max(TRAIL)+1;  % put the new maximum value +1 as an assumed answer
                        else
                            RE1(q) = max(TRAIL); % put the new maximum value only as an assumed answer
                        end
                        flag(q) = 0;  % Set a flag for the assumed value to 0 %%(an answer is found for this unknown)
                        counto = counto +1; % Increment the iteration counter
        end
    end
    
%Looping over the remaining equations to solve them using the assumptions
    for l = 1:(nodes-1)

        RE1(l) = P(l); % putting the right hand side into the unknown

        for q = streets:-1:(nodes)
            RE1(l) = RE1(l) - RE1(q).*x(l,q);  %decrement the values of the assumptions from the right hand side
        end
    
   
        if(abs(RE1(l)) > Constraints(l)) %if the answer exceeded the constraint
            for m = 1:streets %%reset all the values to not found
                flag(m) = 1;
            end 
            break;
        else
            flag(l) = 0; %else, set the unknown to be found
        end
    end
  
  zflag = sum(flag); % summing all the falgs to find out if they are all set to zero

end %%loops back


if (counto < 10000000) %% if counto is less than Million then a first answer is found,
        
        %%%Re-Preform the Previous Steps for finding an answer
        %%% However this time using a new array for the second results RE2
        flag = ones(1, streets);
        zflag = sum(flag);

         for q = streets:-1:(nodes)
            RE2(q) = 1;
         end
        counto = 0;
        while ((zflag ~= 0) && counto < 10000000)

            for q = streets:-1:(nodes)
                tmpo = RE2(q);
                RE2(q) = randi(Constraints(q));

                while ((tmpo > RE2(q) && counto <5000000) || RE2(q) == RE1(q)) 
                    TRAIL = randi(Constraints(q), [50,1]);
                    if(max(TRAIL)+1 <= Constraints(q))
                        RE2(q) = max(TRAIL)+1;
                    else
                        RE2(q) = max(TRAIL);
                    end
                    flag(q) = 0;
                    counto = counto +1;
                end
               if(counto >= 5000000)
                    while ((tmpo < RE2(q)|| RE2(q) == RE1(q)) && counto < 10000000)  %In this loop we consider the case where the new assumption is equal to the first result
                        TRAIL = randi(RE2(q), [50,1]);
                        if(max(TRAIL)+1 <= Constraints(q))
                            RE2(q) = max(TRAIL)+1;
                        else
                            RE2(q) = max(TRAIL);
                        end
                        flag(q) = 0;
                        counto = counto +1;
                    end
               end
            end

            for l = 1:(nodes-1)
                RE2(l) = P(l);

                for q = streets:-1:(nodes)
                    RE2(l) = RE2(l) - RE2(q).*x(l,q);
                end


                if(abs(RE2(l)) > Constraints(l)) 
                    for m = 1:streets
                        flag(m) = 1;
                    end
                    break;
                else
                    flag(l) = 0;
                end
            end

          zflag = sum(flag);

        end


    if(counto < 10000000) %% if the new counto is less than Million then a second answer has been found
        fprintf('The first Result are: ');  %Display the array of the first answer
        disp(RE1);
        fprintf('The Second Result are: '); %Display the array of the second answer
        disp(RE2);
        fprintf('\nFirst Results in order x1, x2, x3:\n');

        %Looping over the first answer array to diplay each value and
        %determine its case: critical/ near critical/ safe
        for h = 1:nodes-1
            disp(RE1(h));
            if(abs(Constraints(h))-abs(RE1(h)) <= 10)
                fprintf('Critical Value\n\n');
            elseif (abs(Constraints(h))-abs(RE1(h)) <= 50)
                fprintf('Near Critical Value\n\n');
            else
                fprintf('Safe Value\n\n');
            end
        end

        fprintf('\nAssumptions are:\n'); %Print out the assumptions

        %Looping over the first assumptions to diplay each value and
        %determine its case: critical/ near critical/ safe
        for h = nodes:streets 
           disp(RE1(h));
            if(abs(Constraints(h))-abs(RE1(h)) <= 10)
                fprintf('Critical Value\n\n');
            elseif (abs(Constraints(h))-abs(RE1(h)) <= 50)
                fprintf('Near Critical Value\n\n');
            else
                fprintf('Safe Value\n\n');
            end
        end

        %%%Repeat the same process to output the second answers
        fprintf('\nSecond Results in order x1, x2, x3...etc:\n');

        for h = 1:nodes-1
            disp(RE2(h))
            if(abs(Constraints(h))-abs(RE2(h)) <= 10)
                fprintf('Critical Value\n\n');
            elseif (abs(Constraints(h))-abs(RE2(h)) <= 50)
                fprintf('Near Critical Value\n\n');
            else
                fprintf('Safe Value\n\n');
            end
        end

        fprintf('\nAssumptions are:\n');

        for h = nodes:streets
           disp(RE2(h))
            if(abs(Constraints(h))-abs(RE2(h)) <= 10)
                fprintf('Critical Value\n\n');
            elseif (abs(Constraints(h))-abs(RE2(h)) <= 50)
                fprintf('Near Critical Value\n\n');
            else
                fprintf('Safe Value\n\n');
            end
        end
    else
        %%%If Only one answer was found
        fprintf('found only one answer within the constraints\nThe Results are: '); %inform the user
        disp(RE1);  % diplay the array of the answer
              
        fprintf('\nResults in order x1, x2, x3:\n');
        
        
        %Looping over the answer array to diplay each value and
        %determine its case: critical/ near critical/ safe
        for h = 1:nodes-1
            disp(RE1(h));
            if(abs(Constraints(h))-abs(RE1(h)) <= (0.05 * (Constraints(h)))
                fprintf('Critical Value\n\n');
            elseif (abs(Constraints(h))-abs(RE1(h)) <= (0.15 * (Constraints(h)))
                fprintf('Near Critical Value\n\n');
            else
                fprintf('Safe Value\n\n');
            end
        end

        fprintf('\nAssumptions are:\n');

        %Looping over the assumptions to diplay each value and
        %determine its case: critical/ near critical/ safe
        for h = nodes:streets
           disp(RE1(h));
            if(abs(Constraints(h))-abs(RE1(h)) <= (0.05 * (Constraints(h))))
                fprintf('Critical Value\n\n');
            elseif (abs(Constraints(h))-abs(RE1(h)) <= (0.15 *(Constraints(h))))
                fprintf('Near Critical Value\n\n');
            else
                fprintf('Safe Value\n\n');
            end
        end
    end
else
    %%%If NO answer was found, inform the user
    fprintf('Cannot find an answer within the given constraints\n');
end
