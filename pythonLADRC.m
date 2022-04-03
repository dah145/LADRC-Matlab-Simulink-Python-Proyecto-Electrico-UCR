function u = pythonLADRC(r, y, param)

    persistent pyLADRC
    
    % Inicializar
    if isempty(pyLADRC)
        pyLADRC = py.LADRC_en_python.LADRC(param.nx, param.b, param.wc, param.wo, param.zo);           
    end    
    
    u = pyLADRC.SalidaControl(r, y);
end

    