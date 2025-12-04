"""
bioconfigme.py - Configuration loader for gene-based analysis pipeline

This module provides functions to load and access configuration values from YAML files.
All Snakemake rules should import this module using:
    sys.path.append("utils")
    from bioconfigme import get_results_dir, get_software_module, get_analysis_value

API Functions:
    - get_results_dir() -> str
    - get_software_module(tool: str) -> str
    - get_software_command(tool: str) -> str
    - get_software_param(tool: str, param: str) -> any
    - get_analysis_value(path: str) -> any
    - get_config_value(path: str) -> any
    - get_dataset(name: str) -> dict
    - get_ref_panel(population: str) -> str
    - get_database(name: str) -> str
"""

import yaml
import os
from pathlib import Path
from typing import Any, Dict, Optional


class ConfigLoader:
    """Singleton configuration loader"""
    _instance = None
    _configs = {}
    
    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(ConfigLoader, cls).__new__(cls)
        return cls._instance
    
    def load_config(self, config_path: str, config_type: str):
        """Load a YAML configuration file"""
        if config_type not in self._configs:
            with open(config_path, 'r') as f:
                self._configs[config_type] = yaml.safe_load(f)
        return self._configs[config_type]
    
    def get(self, config_type: str, path: str = None):
        """Get configuration value by dot-separated path"""
        if config_type not in self._configs:
            raise ValueError(f"Configuration type '{config_type}' not loaded")
        
        config = self._configs[config_type]
        if path is None:
            return config
        
        keys = path.split('.')
        value = config
        for key in keys:
            if isinstance(value, dict):
                value = value.get(key)
            else:
                return None
        return value


# Global instance
_loader = ConfigLoader()


def init_configs(software_path: str, analysis_path: str, config_path: str):
    """Initialize configuration files"""
    _loader.load_config(software_path, 'software')
    _loader.load_config(analysis_path, 'analysis')
    _loader.load_config(config_path, 'config')


def get_results_dir() -> str:
    """Get the results directory path from analysis.yml"""
    return _loader.get('analysis', 'results_dir')


def get_software_module(tool: str) -> str:
    """Get the module name for a software tool"""
    return _loader.get('software', f'{tool}.module')


def get_software_path(tool: str) -> Optional[str]:
    """Get the executable path for a software tool"""
    return _loader.get('software', f'{tool}.path')


def get_software_command(tool: str) -> str:
    """
    Get the command name for a software tool.
    Returns the path if provided, otherwise returns the command name or tool name.
    """
    path = _loader.get('software', f'{tool}.path')
    if path:
        return path
    cmd = _loader.get('software', f'{tool}.command')
    return cmd if cmd else tool


def get_software_param(tool: str, param: str) -> Any:
    """Get a parameter value for a software tool"""
    return _loader.get('software', f'{tool}.params.{param}')


def get_analysis_value(path: str) -> Any:
    """Get a value from analysis.yml by dot-separated path"""
    return _loader.get('analysis', path)


def get_config_value(path: str) -> Any:
    """Get a value from config.yml by dot-separated path"""
    return _loader.get('config', path)


def get_dataset(name: str) -> Optional[Dict]:
    """Get dataset configuration by name"""
    datasets = _loader.get('analysis', 'datasets')
    if datasets:
        for ds in datasets:
            if ds.get('name') == name:
                return ds
    return None


def get_ref_panel(population: str) -> str:
    """Get reference panel path for a population"""
    return _loader.get('analysis', f'ref_panel.{population}')


def get_database(name: str) -> str:
    """Get gene-set database path by name"""
    return _loader.get('analysis', f'databases.{name}')


def get_target_analysis(name: str) -> Optional[Dict]:
    """Get target analysis configuration by name"""
    return _loader.get('analysis', f'target_analysis.{name}')
