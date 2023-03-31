function [solution, Info] = homotopy_solve(solver, args, sEnd)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
kappa_s_times = 0.2;
kappa_s_exp = 1.5;
homotopy_counter = 1;
iterNum = 0;
totalTime_Start = tic;
while true
    % Homotopy (outer) iteration
    HomotopyTime_Start = tic;
    solution = solver('x0', args.x0, 'p',args.p,...
        'lbx', args.lbx, 'ubx', args.ubx,...
        'lbg', args.lbg, 'ubg', args.ubg);
    HomotopyTime = toc(HomotopyTime_Start);
    % check IPOPT status
    if strcmp(solver.stats.return_status, 'Solve_Succeeded')
        % pring information (solver.stats)
        iterNum = iterNum + solver.stats.iter_count;
        msg = ['Homotopy Iter: ', num2str(homotopy_counter), '; ',...
            's: ', num2str(args.p,'%10.2e'), '; ',...
            'Cost: ', num2str(full(solution.f),'%10.2e'), '; ',...
            'Ipopt IterNum: ', num2str(solver.stats.iter_count), '; ',...
            'Time: ', num2str(1000*HomotopyTime,'%10.2e'), ' ms'];
        disp(msg)
    else
        % IPOPT fails
        disp('IPOPT fails')
        homotopy_status = 0;
        break
    end
    % check termination of homotopy
    if (args.p == sEnd) && (strcmp(solver.stats.return_status, 'Solve_Succeeded'))
        % success
        homotopy_status = 1;
        break
    else
        % update x0 and s for next homotopy iteration
        args.x0 = solution.x;
        s_trail = min([kappa_s_times .* args.p, args.p.^kappa_s_exp]);
        args.p = max([s_trail, sEnd]);
        homotopy_counter = homotopy_counter + 1;        
    end
end
totalTime = toc(totalTime_Start);

Info = struct('iterNum', iterNum, 'totalTime', totalTime, 'homotopy_status', homotopy_status);

end

