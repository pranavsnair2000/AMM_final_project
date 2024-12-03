
function coppelia(init_pos, Thetas ,Tspan)
    SimStarted = false;
    % disp(size(Thetas))
    % disp(size(Tspan))

    % %% Connect to Coppeliasim
    client = RemoteAPIClient();         % Create a client object
    sim = client.getObject('sim');      % Reference the Sim Object
    
    defaultIdleFps = sim.getInt32Param(sim.intparam_idle_fps); % Get default parameter for fps
    sim.setInt32Param(sim.intparam_idle_fps, 0); % Control the fps
    
    % world = sim.getObject('./Floor');    
    
    % get object handles
    j1= sim.getObject('./j1');
    j2= sim.getObject('./j2');
    j3= sim.getObject('./j3');
    j4= sim.getObject('./j4');
    j5= sim.getObject('./j5');

    
    
    
    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % client.setStepping(true); 
    sim.startSimulation(); % Start Simulation
    SimStarted = true;
    t = sim.getSimulationTime(); % Keep track of time
    
    % initial position
    sim.setJointTargetPosition(j1, init_pos(1));
    sim.setJointTargetPosition(j2, init_pos(2));
    sim.setJointTargetPosition(j3, init_pos(3));
    sim.setJointTargetPosition(j4, init_pos(4));
    sim.setJointTargetPosition(j5, init_pos(5));
    
    pause(1)
    
    for i = 1:length(Tspan)
        sim.setJointTargetPosition(j1, Thetas(i,1));
        sim.setJointTargetPosition(j2, Thetas(i,2));
        sim.setJointTargetPosition(j3, Thetas(i,3));
        sim.setJointTargetPosition(j4, Thetas(i,4));
        sim.setJointTargetPosition(j5, Thetas(i,5));
        
        % pause(20/length(Tspan))
    end

    pause(5)

    % Destructor
    
    if SimStarted
        sim.stopSimulation();    
    end
    
    sim.setInt32Param(sim.intparam_idle_fps, defaultIdleFps);
end