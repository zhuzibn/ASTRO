function fields = astro_evaluate_field_terms(mmxtmp, mmytmp, mmztmp, ...
        AsimnextL, AsimpreviousL, AsimnextW, AsimpreviousW, muigpu, ...
        Ksim, Dsim, Hext, bc)
%ASTRO_EVALUATE_FIELD_TERMS Evaluate deterministic effective-field terms.
%   Returns exchange, anisotropy, DMI, external, and total fields in tesla.
%   Thermal and dipole terms stay outside this helper because they are not
%   deterministic field terms in the production script.

mmxnextL = circshift(mmxtmp, [0, -1]);%length direction
mmynextL = circshift(mmytmp, [0, -1]);
mmznextL = circshift(mmztmp, [0, -1]);
mmxpreviousL = circshift(mmxtmp, [0, 1]);
mmypreviousL = circshift(mmytmp, [0, 1]);
mmzpreviousL = circshift(mmztmp, [0, 1]);

mmxnextW = circshift(mmxtmp, [1, 0]);%width direction
mmynextW = circshift(mmytmp, [1, 0]);
mmznextW = circshift(mmztmp, [1, 0]);
mmxpreviousW = circshift(mmxtmp, [-1, 0]);
mmypreviousW = circshift(mmytmp, [-1, 0]);
mmzpreviousW = circshift(mmztmp, [-1, 0]);

if bc%not periodic condition
    mmxnextL(:, end) = 0;mmynextL(:, end) = 0;mmznextL(:, end) = 0;
    mmxpreviousL(:, 1) = 0;mmypreviousL(:, 1) = 0;mmzpreviousL(:, 1) = 0;
    mmxnextW(1, :) = 0;mmynextW(1, :) = 0;mmznextW(1, :) = 0;
    mmxpreviousW(end, :) = 0;mmypreviousW(end, :) = 0;mmzpreviousW(end, :) = 0;
else%periodic condition
    %do nothing
end

fields.exchange.x = -(AsimnextL.*mmxnextL + AsimpreviousL.*mmxpreviousL + ...
    AsimnextW.*mmxnextW + AsimpreviousW.*mmxpreviousW)./muigpu;%[T]
fields.exchange.y = -(AsimnextL.*mmynextL + AsimpreviousL.*mmypreviousL + ...
    AsimnextW.*mmynextW + AsimpreviousW.*mmypreviousW)./muigpu;
fields.exchange.z = -(AsimnextL.*mmznextL + AsimpreviousL.*mmzpreviousL + ...
    AsimnextW.*mmznextW + AsimpreviousW.*mmzpreviousW)./muigpu;

fields.anisotropy.x = zeros(size(fields.exchange.x), 'like', fields.exchange.x);%anisotropy
fields.anisotropy.y = zeros(size(fields.exchange.x), 'like', fields.exchange.x);
fields.anisotropy.z = 2*Ksim./muigpu.*mmztmp;%[T]

fields.dmi.x = Dsim./muigpu.*(-mmznextL + mmzpreviousL);%[T]
fields.dmi.y = Dsim./muigpu.*(-mmznextW + mmzpreviousW);
fields.dmi.z = Dsim./muigpu.*(mmxnextL - mmxpreviousL + mmynextW - mmypreviousW);

fields.external.x = Hext(1) + zeros(size(fields.exchange.x), 'like', fields.exchange.x);
fields.external.y = Hext(2) + zeros(size(fields.exchange.x), 'like', fields.exchange.x);
fields.external.z = Hext(3) + zeros(size(fields.exchange.x), 'like', fields.exchange.x);

fields.total.x = fields.exchange.x + fields.anisotropy.x + fields.dmi.x + fields.external.x;
fields.total.y = fields.exchange.y + fields.anisotropy.y + fields.dmi.y + fields.external.y;
fields.total.z = fields.exchange.z + fields.anisotropy.z + fields.dmi.z + fields.external.z;
end
