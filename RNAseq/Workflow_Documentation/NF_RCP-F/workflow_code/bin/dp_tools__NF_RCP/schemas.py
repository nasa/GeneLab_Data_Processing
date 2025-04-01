import pandas as pd

class runsheet: # Bad casing since we will use the class definition itself for all static methods
    
    @staticmethod
    def check_single_value(column: pd.Series, error_msg: str, errors: list[str]) -> None:
        if len(column.unique()) != 1:
            errors.append(error_msg)
    
    @staticmethod
    def check_read2_path_populated_if_paired_end(df: pd.DataFrame, errors: list[str]) -> None:
        if (("read2_path" in df.columns and df['paired_end'][0] == True) or
           ("read2_path" not in df.columns and df['paired_end'][0] == False)):
            return
        else:
            errors.append("Expected 'read2_path' to be populated only if paired_end is True")
    
    @staticmethod
    def validate(df_runsheet: pd.DataFrame) -> bool:
        errors = []
        
        # Check for single value in specified columns
        runsheet.check_single_value(df_runsheet['has_ERCC'], "Dataset level columns do NOT contain one unique value for 'has_ERCC'", errors)
        runsheet.check_single_value(df_runsheet['organism'], "Dataset level columns do NOT contain one unique value for 'organism'", errors)
        runsheet.check_single_value(df_runsheet['paired_end'], "Dataset level columns do NOT contain one unique value for 'paired_end'", errors)
        
        # Check for 'read2_path' population if paired_end is True
        runsheet.check_read2_path_populated_if_paired_end(df_runsheet, errors)
        
        if errors:
            raise ValueError("\n".join(errors))
        else:
            return True