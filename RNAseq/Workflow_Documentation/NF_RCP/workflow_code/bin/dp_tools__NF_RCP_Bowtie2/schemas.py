import pandas as pd

class runsheet: # Bad casing since we will use the class definition itself for all static methods
    
    @staticmethod
    def check_single_value(column: pd.Series, error_msg: str, errors: list[str]) -> None:
        if len(column.unique()) != 1:
            errors.append(error_msg)
    
    @staticmethod
    def check_read2_path_populated_if_paired_end(df: pd.DataFrame, errors: list[str]) -> None:
        if (("read2_path" in df.columns and df['paired_end'].iloc[0] == True) or
            ("read2_path" not in df.columns and df['paired_end'].iloc[0] == False)):
            return
        else:
            errors.append("Expected 'read2_path' to be populated only if paired_end is True")
    
    @staticmethod
    def coerce_boolean_columns(df: pd.DataFrame) -> pd.DataFrame:
        """Coerce `Has Tech Reps` column to boolean type"""
        if 'Has Tech Reps' in df.columns:
            # Handle string boolean values
            df['Has Tech Reps'] = df['Has Tech Reps'].map(
                lambda x: str(x).lower() == 'true' if pd.notna(x) else False
            )
        
        return df
    
    @staticmethod
    def check_source_name_and_has_tech_reps_dependency(df: pd.DataFrame, errors: list[str]) -> None:
        if "Has Tech Reps" in df.columns and "Source Name" not in df.columns:
            errors.append("'Source Name' column must be present when 'Has Tech Reps' column is used")

    @staticmethod
    def validate(df_runsheet: pd.DataFrame) -> bool:
        errors = []
        
        # Coerce boolean columns first
        df_runsheet = runsheet.coerce_boolean_columns(df_runsheet)
        
        # Check for single value in specified columns
        runsheet.check_single_value(df_runsheet['has_ERCC'], "Dataset level columns do NOT contain one unique value for 'has_ERCC'", errors)
        runsheet.check_single_value(df_runsheet['organism'], "Dataset level columns do NOT contain one unique value for 'organism'", errors)
        runsheet.check_single_value(df_runsheet['paired_end'], "Dataset level columns do NOT contain one unique value for 'paired_end'", errors)
        
        # Check for 'read2_path' population if paired_end is True
        runsheet.check_read2_path_populated_if_paired_end(df_runsheet, errors)
        
        # Check for Source Name and Has Tech Reps dependency
        runsheet.check_source_name_and_has_tech_reps_dependency(df_runsheet, errors)
        
        if errors:
            raise ValueError("\n".join(errors))
        else:
            return True