function Penalty = PsiPenalty(x_cur,y_cur)
    E = 0;
    penetration_x = max(abs(x_cur - 1) - 1,0);
    penetration_y = max(abs(y_cur - 1/2) - 1/2,0);
    penetration2 = penetration_x^2 + penetration_y^2;
    Penalty = penetration2 * E * 1/2;
end
