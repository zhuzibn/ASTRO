function atomtype_ = astro_validate_atom_distribution( ...
        atomtype_, natomW, natomL, compositionn)
%ASTRO_VALIDATE_ATOM_DISTRIBUTION Validate deterministic RE/TM site matrices.

if ~isequal(size(atomtype_), [natomW, natomL])
    error('ASTRO:InvalidAtomDistribution', ...
        'atomtype_ must have size [%d, %d].', natomW, natomL);
end

if ~all(atomtype_(:) == 0 | atomtype_(:) == 1)
    error('ASTRO:InvalidAtomDistribution', ...
        'atomtype_ must contain only 0 (TM/FeCo) and 1 (RE/Gd).');
end

if nargin >= 4 && ~isempty(compositionn)
    expected_gd_count = round(natomW*natomL*compositionn);
    actual_gd_count = nnz(atomtype_ == 1);
    if actual_gd_count ~= expected_gd_count
        error('ASTRO:InvalidAtomDistribution', ...
            ['atomtype_ contains %d RE/Gd sites; expected %d for ' ...
            'compositionn=%g.'], actual_gd_count, expected_gd_count, ...
            compositionn);
    end
end
end
