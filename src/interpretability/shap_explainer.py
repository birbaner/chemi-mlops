import sys
from pathlib import Path
import numpy as np
import shap
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestRegressor

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))


FEATURE_NAMES = [
    "MolWt",
    "LogP", 
    "TPSA",
    "NumHDonors",
    "NumHAcceptors",
    "NumRotatableBonds",
    "RingCount"
]


class ModelInterpreter:
    """Provides model interpretability via SHAP and uncertainty quantification."""
    
    def __init__(self, model: RandomForestRegressor, X_train: np.ndarray):
        """
        Initialize the interpreter with a trained model and training data.
        
        Args:
            model: Trained RandomForestRegressor
            X_train: Training features used to fit the explainer
        """
        self.model = model
        self.explainer = shap.TreeExplainer(model)
        # Use a subset of training data for faster SHAP computation
        self.X_background = shap.sample(X_train, min(100, len(X_train)))
    
    def get_feature_importance(self) -> dict:
        """Get global feature importance from the random forest model."""
        importance = self.model.feature_importances_
        return {
            name: float(imp) 
            for name, imp in zip(FEATURE_NAMES, importance)
        }
    
    def explain_prediction(self, X: np.ndarray):
        """
        Explain a single prediction using SHAP values.
        
        Args:
            X: Feature vector (1, n_features)
            
        Returns:
            dict with SHAP values and base value
        """
        shap_values = self.explainer.shap_values(X)
        return {
            "shap_values": {
                name: float(val) 
                for name, val in zip(FEATURE_NAMES, shap_values[0])
            },
            "base_value": float(self.explainer.expected_value),
            "prediction": float(self.model.predict(X)[0])
        }
    
    def get_prediction_interval(self, X: np.ndarray, confidence=0.95) -> tuple:
        """
        Estimate prediction interval using individual tree predictions.
        
        Args:
            X: Feature vector (1, n_features)
            confidence: Confidence level (default 0.95)
            
        Returns:
            (lower_bound, upper_bound, std_dev)
        """
        # Get predictions from all trees
        tree_predictions = np.array([
            tree.predict(X)[0] 
            for tree in self.model.estimators_
        ])
        
        mean_pred = tree_predictions.mean()
        std_pred = tree_predictions.std()
        
        # Calculate confidence interval
        alpha = 1 - confidence
        lower = np.percentile(tree_predictions, alpha/2 * 100)
        upper = np.percentile(tree_predictions, (1 - alpha/2) * 100)
        
        return float(lower), float(upper), float(std_pred)
    
    def plot_shap_waterfall(self, X: np.ndarray, feature_values: dict) -> plt.Figure:
        """
        Create a SHAP waterfall plot showing feature contributions.
        
        Args:
            X: Feature vector (1, n_features)
            feature_values: Dict mapping feature names to values
            
        Returns:
            matplotlib Figure
        """
        shap_values = self.explainer.shap_values(X)
        base_value = self.explainer.expected_value
        
        # Create waterfall plot
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Sort by absolute SHAP value
        indices = np.argsort(np.abs(shap_values[0]))[::-1]
        
        y_pos = np.arange(len(FEATURE_NAMES))
        colors = ['#ff0051' if val > 0 else '#008bfb' for val in shap_values[0][indices]]
        
        ax.barh(y_pos, shap_values[0][indices], color=colors)
        ax.set_yticks(y_pos)
        ax.set_yticklabels([f"{FEATURE_NAMES[i]} = {X[0][i]:.2f}" for i in indices])
        ax.set_xlabel('SHAP value (impact on prediction)')
        ax.set_title(f'Feature Contributions to Prediction\nBase value: {base_value:.3f}')
        ax.axvline(x=0, color='black', linestyle='-', linewidth=0.5)
        
        plt.tight_layout()
        return fig
