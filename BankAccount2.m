classdef BankAccount
    properties (Access=private)
        Balance
        historyLog
        Pin
     
    end

    methods
        function obj=BankAccount(initialBalance,pin)
            if initialBalance>0
                obj.Balance=initialBalance;
                obj.Pin=pin;
                obj.historyLog={};
            else
                error('Initial balance must be non-negative.')
            end
        end

        function obj=deposit(obj,amount)
            if enterPin(obj)
            if amount>0 && obj.Balance+amount<10000
                obj.Balance=obj.Balance+amount;
                obj.historyLog=[obj.historyLog; sprintf('Deposited $%.2f',amount)];
            elseif obj.Balance+amount>=10000
                error('Balance cannot be higher than $10000.')
            else
                error('Deposit amount has to be positive.')
            end
            else
                error('Invalid PIN');
            end
        end

        function [obj, amountWithdrawn]=withdraw(obj,amount)
            if enterPin(obj)
            if amount<=obj.Balance && amount>0 && amount<1000
                obj.Balance=obj.Balance-amount;
                amountWithdrawn=amount;
                obj.historyLog=[obj.historyLog; sprintf('Withdrawn $%.2f',amount)];
            elseif amount>=1000
                error('Cannot withdraw more than $1000 at once.')
            else
                error('Insufficient funds or invalid amount.')
            end
            else
                error('Invalid PIN');
            end
        end

        function getBalance(obj)
            if enterPin(obj)
             disp(obj.Balance);
            else
                error('Invalid PIN')
            end

        end

        function [obj1,obj2]=transferMoney(obj1,obj2,amount)
           if enterPin(obj1) 
            if obj1.Balance<amount && amount<10
                error('Insufficient funds or amount to be transferred is less than $10.')
            elseif amount>1000
                error('Amount cannot be higher than $1000.')
            else
                obj1.Balance=obj1.Balance-amount;
                obj2.Balance=obj2.Balance+amount;
                obj1.historyLog=[obj1.historyLog; sprintf('Transferred $%.2f',amount)];
                obj2.historyLog=[obj2.historyLog; sprintf('Received $%.2f',amount)];
            end
           else
               error('Invalid PIN');
           end
        end
        function ok=enterPin(obj)
            prompt='Enter PIN: ';
            x=input(prompt);
            ok=(x==obj.Pin);
        end
        function getHistory(obj)
            if enterPin(obj)
            disp(obj.historyLog);
            else
                error('Invalid PIN');
            end
        end
    end
end