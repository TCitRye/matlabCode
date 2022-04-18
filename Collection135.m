function [lifeValue, n] = Collection135(ATT, ATTBouns, Hp1, MGDefence, MGBonus)
    n = 0;
    Hp(1) = Hp1;
    Damage = ATT*(1+ATTBonus)*(1+MGBonus)*(1-MGDefence);

    while Hp >= 0
        Hp(n+2) = Hp(n+1) - 0.03*Hp(n+1) - Damage;
        n=n+1;
    end

    n = [0:1:n];

    figure
    plot(n, Hp)
    xlabel('n')
    ylabel('Hp')
