####################################################################################################
# # Copyright (C) 2023-Present - Daniel Charytonowicz - All Rights Reserved
# # Contact: daniel.charytonowicz@icahn.mssm.edu
# ###################################################################################################

# State watcher: Watch for task submission and run state
from typing import Optional

import time
import progressbar
from functools import partial

from ._api import get_task_state, get_run_exists
from ._api import get_run_state, get_run_progress, get_run_result
from ._apiutils import waitfor

def wait_for_submission(
    run_id : str,
    runtype : str,
    token : Optional[str] = None,
    poll_interval : int = 1,
    timeout : int = 3600,
    showprogress : bool = True,
    ) -> None:
    """\
    
    Wait for Submission
    
    Polls unicell API for queue status updates at regular intervals
    and waits until run is actually submitted.
    
    Params
    -------
    run_id
        The run_id
    runtype
        The runtype submitted, used to poll the correct
        task queue to check for run submisssion.
    token
        Optional token, inferred if not given
    poll_interval
        Interval to poll queue for submission, default 4 seconds
    timeout
        Time in seconds to wait before returning exception.
    showprogress
        Show progressbar
        
    Returns
    -------
    
    """
    
    # Get initial state of task
    state_task = get_task_state(run_id, runtype, token)
    
    # Check if run exists yet
    run_exists = get_run_exists(run_id, token)['exists']
    
    # Get initial time
    start_time = time.time()
    
    # Get initial elapsed time
    elapsed_time = time.time() - start_time
    
    # Get counter for progressbar
    counter = 0
    
    # Create prefix definition
    prefix = "Waiting For Submission : {variables.task_state} | Queue Size : {variables.queue_size} | "
    
    # Only show if requested
    if showprogress:

        pbar = progressbar.ProgressBar(
                variables = {"task_state" : state_task['task_state'],
                             "queue_size" : state_task['queue_size']},
                prefix = prefix,
                max_value = None, initial_value = counter)
        pbar.update(counter, task_state = state_task['task_state'],
                   queue_size = state_task['queue_size'])
        
    # Loop while task is not yet running
    while (not (state_task["task_state"] == "RUNNING")) and (not run_exists):

        # Wait interval
        time.sleep(poll_interval)

        # Update task state and run status
        state_task = get_task_state(run_id, runtype, token)
        run_exists = get_run_exists(run_id, token)['exists']
        
        # Increment counter
        counter += 1
        
        # Only update if requested
        if showprogress:
            
            # Update progress bar
            pbar.update(counter, 
                task_state = state_task['task_state'],
                queue_size = state_task['queue_size'])
        
        # Update elapsed time
        elapsed_time = time.time() - start_time
        
        # Check for timeout
        if elapsed_time > timeout:
            raise Exception("Task state timeout, please retry.")
        
    
    # Finish progressbar if present
    if showprogress: 
        pbar.update(task_state = "RUNNING")
        pbar.finish()
        

def wait_for_completion(
    run_id : str,
    token : Optional[str] = None,
    poll_interval : int = 1,
    timeout : int = 3600,
    showprogress : bool = True,
    pbar_key : str = "progress",
    desc : Optional[str] = None,
    ) -> None:
    """\
    Wait for Completion
    
    Polls unicell API for status updates at regular intervals about given run until
    run reports back complete or timeout reached.
    First waits for task to be submitted, then checks run progress. Either state is
    completed or killed, either way we can then end the task
    
    Params
    -------
    run_id
        The run_id
    token
        Optional token, inferred if not given
    poll_interval
        Interval to poll, default 2 seconds.
    timeout
        Time in seconds to wait before forcing return on run
    showprogress
        Get progressbar from backend
    pbar_key
        Key for progressbar on backend
    desc
        Description for progressbar, default None
    Returns
    -------
    
    """
    
    # Get initial state
    state = get_run_state(run_id, token)
    
    # Get initial time
    start_time = time.time()
    
    # Get initial elapsed time
    elapsed_time = time.time() - start_time

    # Create a progressbar if desired
    if showprogress:
        
        # Get initial progress state using a waitfor wrapper
        runprogress = waitfor(partial(get_run_progress, run_id = run_id, token = token))
        
        pbar = progressbar.ProgressBar(
            prefix = f"{desc} | " if desc else None,
            max_value = runprogress['total_steps'],
            initial_value = runprogress['current_step'])
        pbar.update(runprogress['current_step'])
    else:
        pbar = None
        

    # Loop while state is not completed or timeout not reached
    while (not state['completed']) and (elapsed_time < timeout) and (not state['kill']):

        # Wait interval
        time.sleep(poll_interval)

        # Update state and progress
        state = get_run_state(run_id, token)

        # Update progress if desired
        if showprogress:
            
            # We wrap in a waitfor to avoid race condtion where an error gets triggered
            # during a write from the backend which causes the run read to be null
            runprogress = waitfor(partial(get_run_progress, run_id = run_id, token = token))
            pbar.update(runprogress['current_step'])
    
    # Finish progressbar if present
    if showprogress: 
        pbar.finish()
        
    # Once complete get run result if run was not killed
    if not state['kill']:
        return get_run_result(run_id, token)