from flask_sqlalchemy import SQLAlchemy
from flask_login import UserMixin
from datetime import datetime
import json

db = SQLAlchemy()

class User(UserMixin, db.Model):
    """
    User model for authentication and profile management.
    Prepared for future expansion with OAuth/ORCID.
    """
    __tablename__ = 'users'
    
    id = db.Column(db.Integer, primary_key=True)
    username = db.Column(db.String(80), unique=True, nullable=False)
    email = db.Column(db.String(120), unique=True, nullable=False)
    password_hash = db.Column(db.String(128)) # To be implemented with werkzeug.security
    created_at = db.Column(db.DateTime, default=datetime.utcnow)
    
    workspaces = db.relationship('Workspace', backref='owner', lazy=True)

class Workspace(db.Model):
    """
    Individual research workspaces to group analysis results.
    """
    __tablename__ = 'workspaces'
    
    id = db.Column(db.Integer, primary_key=True)
    title = db.Column(db.String(200), default="Untitled Analysis")
    user_id = db.Column(db.Integer, db.ForeignKey('users.id'), nullable=False)
    created_at = db.Column(db.DateTime, default=datetime.utcnow)
    updated_at = db.Column(db.DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)
    
    results = db.relationship('AnalysisResult', backref='workspace', lazy=True)

class AnalysisResult(db.Model):
    """
    Cached results for expensive bioinformatics computations.
    """
    __tablename__ = 'analysis_results'
    
    id = db.Column(db.Integer, primary_key=True)
    tool_type = db.Column(db.String(50), nullable=False) # e.g., 'alignment', 'primer'
    input_params = db.Column(db.Text) # JSON string of inputs
    result_data = db.Column(db.Text)  # JSON string of outputs
    workspace_id = db.Column(db.Integer, db.ForeignKey('workspaces.id'), nullable=True)
    created_at = db.Column(db.DateTime, default=datetime.utcnow)

    def set_params(self, data):
        self.input_params = json.dumps(data)

    def get_params(self):
        return json.loads(self.input_params) if self.input_params else {}

    def set_results(self, data):
        self.result_data = json.dumps(data)

    def get_results(self):
        return json.loads(self.result_data) if self.result_data else {}
