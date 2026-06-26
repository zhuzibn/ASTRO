function [AsimnextL, AsimpreviousL, AsimnextW, AsimpreviousW] = ...
        astro_exchange_neighbor_coefficients(atomtype_, bc, Jgdgd, Jfefe, Jfegd)
%ASTRO_EXCHANGE_NEIGHBOR_COEFFICIENTS Build nearest-neighbor exchange links.
%   atomtype_ == 1 is RE/Gd and atomtype_ == 0 is TM/FeCo.  The first matrix
%   dimension is physical width W and the second is physical length L.

[natomW, natomL] = size(atomtype_);

AsimnextL = zeros(natomW, natomL);
AsimpreviousL = zeros(natomW, natomL);
AsimnextW = zeros(natomW, natomL);
AsimpreviousW = zeros(natomW, natomL);

for ctL = 1:natomL
    for ctW = 1:natomW
        if atomtype_(ctW, ctL) == 1%RE
            if ctL == natomL
                if bc
                    AsimnextL(ctW, ctL) = 0;
                else
                    AsimnextL(ctW, ctL) = Jgdgd*(atomtype_(ctW, ctL) == atomtype_(ctW, 1)) + ...
                        Jfegd*(~atomtype_(ctW, ctL) == atomtype_(ctW, 1));
                end
            else
                if atomtype_(ctW, ctL + 1) == 1%the other atom is RE
                    AsimnextL(ctW, ctL) = Jgdgd;
                else%the other atom is TM
                    AsimnextL(ctW, ctL) = Jfegd;
                end
            end

            if ctL == 1
                if bc
                    AsimpreviousL(ctW, ctL) = 0;
                else
                    AsimpreviousL(ctW, ctL) = Jgdgd*(atomtype_(ctW, ctL) == atomtype_(ctW, natomL)) + ...
                        Jfegd*(~atomtype_(ctW, ctL) == atomtype_(ctW, natomL));
                end
            else
                if atomtype_(ctW, ctL - 1) == 1%the other atom is RE
                    AsimpreviousL(ctW, ctL) = Jgdgd;
                else%the other atom is TM
                    AsimpreviousL(ctW, ctL) = Jfegd;
                end
            end

            if ctW == natomW
                if bc
                    AsimpreviousW(ctW, ctL) = 0;
                else
                    AsimpreviousW(ctW, ctL) = Jgdgd*(atomtype_(ctW, ctL) == atomtype_(1, ctL)) + ...
                        Jfegd*(~atomtype_(ctW, ctL) == atomtype_(1, ctL));
                end
            else
                if atomtype_(ctW + 1, ctL) == 1%the other atom is RE
                    AsimpreviousW(ctW, ctL) = Jgdgd;
                else%the other atom is TM
                    AsimpreviousW(ctW, ctL) = Jfegd;
                end
            end

            if ctW == 1
                if bc
                    AsimnextW(ctW, ctL) = 0;
                else
                    AsimnextW(ctW, ctL) = Jgdgd*(atomtype_(ctW, ctL) == atomtype_(natomW, ctL)) + ...
                        Jfegd*(~atomtype_(ctW, ctL) == atomtype_(natomW, ctL));
                end
            else
                if atomtype_(ctW - 1, ctL) == 1%the other atom is RE
                    AsimnextW(ctW, ctL) = Jgdgd;
                else%the other atom is TM
                    AsimnextW(ctW, ctL) = Jfegd;
                end
            end

        else%local atom is TM
            if ctL == natomL
                if bc
                    AsimnextL(ctW, ctL) = 0;
                else
                    AsimnextL(ctW, ctL) = Jfefe*(atomtype_(ctW, ctL) == atomtype_(ctW, 1)) + ...
                        Jfegd*(~atomtype_(ctW, ctL) == atomtype_(ctW, 1));
                end
            else
                if atomtype_(ctW, ctL + 1) == 1%the other atom is RE
                    AsimnextL(ctW, ctL) = Jfegd;
                else%the other atom is TM
                    AsimnextL(ctW, ctL) = Jfefe;
                end
            end

            if ctL == 1
                if bc
                    AsimpreviousL(ctW, ctL) = 0;
                else
                    AsimpreviousL(ctW, ctL) = Jfefe*(atomtype_(ctW, ctL) == atomtype_(ctW, natomL)) + ...
                        Jfegd*(~atomtype_(ctW, ctL) == atomtype_(ctW, natomL));
                end
            else
                if atomtype_(ctW, ctL - 1) == 1%the other atom is RE
                    AsimpreviousL(ctW, ctL) = Jfegd;
                else%the other atom is TM
                    AsimpreviousL(ctW, ctL) = Jfefe;
                end
            end

            if ctW == natomW
                if bc
                    AsimpreviousW(ctW, ctL) = 0;
                else
                    AsimpreviousW(ctW, ctL) = Jfefe*(atomtype_(ctW, ctL) == atomtype_(1, ctL)) + ...
                        Jfegd*(~atomtype_(ctW, ctL) == atomtype_(1, ctL));
                end
            else
                if atomtype_(ctW + 1, ctL) == 1%the other atom is RE
                    AsimpreviousW(ctW, ctL) = Jfegd;
                else%the other atom is TM
                    AsimpreviousW(ctW, ctL) = Jfefe;
                end
            end

            if ctW == 1
                if bc
                    AsimnextW(ctW, ctL) = 0;
                else
                    AsimnextW(ctW, ctL) = Jfefe*(atomtype_(ctW, ctL) == atomtype_(natomW, ctL)) + ...
                        Jfegd*(~atomtype_(ctW, ctL) == atomtype_(natomW, ctL));
                end
            else
                if atomtype_(ctW - 1, ctL) == 1%the other atom is RE
                    AsimnextW(ctW, ctL) = Jfegd;
                else%the other atom is TM
                    AsimnextW(ctW, ctL) = Jfefe;
                end
            end
        end
    end
end
end
