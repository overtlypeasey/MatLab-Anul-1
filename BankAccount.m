classdef BankAccount
    properties 
        Balance
        AccesLog={}
    end
    methods
        function obj = BankAccount(initialBalance)
            if initialBalance >= 0
                obj.Balance = initialBalance;
            else
                error('Initial balance must be positive.');
            end
        end
        
        
        function getBalance(obj)
        disp(['Account balance is: ', num2str(obj.Balance)])
        end
        
        
        function [obj, amountW] = withdraw(obj, amount)
            if obj.Balance >= amount
                obj.Balance = obj.Balance - amount;
                amountW = amount;
                disp('Withdraw succesful')
            else
                error('Insufficient funds')
                amountW = 0;
            end
            operation = {[datetime('now')], ['withdrawal']};
            obj.AccesLog = {[obj.AccesLog]; [operation]};
        end

        function obj = deposit(obj, amount)
            if amount > 0
                obj.Balance = obj.Balance + amount;
                disp('Deposit succesful')
            else
                error('Amount must be positive')
            end
        end

        function [b1,b2]=transferMoney(b1,b2,amount)
            if amount>0
                [b1,withdrawn] = withdraw(b1, amount);
                if withdrawn>0
                    b2 = deposit(b2, withdrawn);
                    disp(['Transferred ', num2str(withdrawn), ' from account 1 to account 2'])
                else 
                    disp('transfer failed')
                end 
            else 
                disp('transfer amount should be positive')
            end 
        end

        function getAccesLog(obj)
            disp(obj.AccesLog{1, :});
        end
    end
end