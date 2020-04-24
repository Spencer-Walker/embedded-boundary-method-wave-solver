function  iP = species_IP(species)
% returns the ionization potential for various species, based on prior DFT
% calculations
if isnumeric(species)
    iP = species;
else
    switch species
        case 'H'
            iP=0.5;
        case 'He'
            iP=0.918;
        case 'Li+'
            iP=2.793;
        case 'Li'
            iP=.2;
        case 'Li-'
            iP=0.015;
        case 'Be'
            iP=0.308;
        case 'F'
            iP=0.709;
        case 'Ne+'
            iP=1.577;
        case 'F-'
            iP=0.16;
        case 'Ne'
            iP=0.821;
        case 'Na+'
            iP=1.737;
        case 'Mg'
            iP=0.281;
        case 'Cl-'
            iP=0.139;
        case 'Ar'
            iP=0.575;
        case 'Ar+'
            iP=1.04919;
        case 'Ar3+'
            iP=2.13571;
        case 'Ar4+'
            iP=2.75235;
        case 'Ar5+'
            iP=3.42072;
        case 'K+'
            iP=1.153;
        case 'Ca'
            iP=0.234;
        case 'Kr'
            iP=0.523;
    end
end
end