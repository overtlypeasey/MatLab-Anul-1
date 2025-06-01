clear
clc
Account = BankAccount(100);
Account = Account.deposit(50);
[Account, amount] = Account.withdraw(100);
Ac2 = BankAccount(0);
[Account, Ac2] = transferMoney(Account, Ac2, 50);
Account.getBalance();
Account.getAccesLog;