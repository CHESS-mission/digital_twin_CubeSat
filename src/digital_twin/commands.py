"""Simple command system for real-time simulation control."""

from typing import Any, Dict


class CommandProcessor:
    """Processes commands and applies them to the simulation."""
    
    def __init__(self, simulation):
        self.simulation = simulation
    
    def execute_command(self, command: Dict[str, Any]) -> tuple[bool, str]:
        """Execute a command by calling the appropriate simulation method.
        
        Args:
            command: Dictionary with 'command' (str) and 'params' (dict) keys
            
        Returns:
            Tuple of (success: bool, message: str)
        """
        try:
            command_name = command.get("command")
            params = command.get("params", {})
            
            if command_name == "set_mode":
                self.simulation.switch_algo.set_mode(params["mode"])
                return True, f"Mode set to {params['mode']}"
            
            elif command_name == "uplink_safe_mode":
                self.simulation.spacecraft.get_telecom().add_uplink_safe_mode(params)
                return True, f"Uplink safe mode commands added: {params}"
            
            else:
                raise ValueError(f"Unknown command: {command_name}")
                
        except Exception as e:
            error_msg = f"Error: {str(e)}"
            return False, error_msg
