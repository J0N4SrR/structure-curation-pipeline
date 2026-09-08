import React, { useEffect, useMemo } from 'react';
import {
  Streamlit,
  withStreamlitConnection,
  ComponentProps,
} from 'streamlit-component-lib';
import ReactFlow, {
  Background,
  Controls,
  Handle,
  Position,
  MarkerType,
} from 'reactflow';
import 'reactflow/dist/style.css';
import { Beaker, ShieldCheck, Combine, Ban, TestTube, CheckCircle } from 'lucide-react';

// Custom Node for the pipeline
const PipelineNode = ({ data }: any) => {
  // Map python glyphs to lucide icons
  const iconMap: any = {
    "🧪": <TestTube size={18} />,
    "⚙️": <Beaker size={18} />,
    "🛡️": <ShieldCheck size={18} />,
    "🧩": <Combine size={18} />,
    "🚫": <Ban size={18} />,
    "✅": <CheckCircle size={18} />
  };
  
  const Icon = iconMap[data.glyph] || <Beaker size={18} />;

  return (
    <div 
      className="pipeline-node"
      style={{
        padding: '10px 14px',
        borderRadius: '8px',
        background: '#ffffff',
        border: `2px solid ${data.color}`,
        minWidth: '200px',
        boxShadow: '0 4px 6px rgba(0,0,0,0.05)',
        display: 'flex',
        flexDirection: 'column',
        gap: '6px'
      }}
    >
      <Handle type="target" position={Position.Top} style={{ background: '#555' }} />
      
      <div style={{ display: 'flex', alignItems: 'center', gap: '8px' }}>
        <div style={{ color: data.color }}>{Icon}</div>
        <strong style={{ fontSize: '13px', color: '#1f2937' }}>{data.label}</strong>
      </div>
      
      <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-end', marginTop: '4px' }}>
        <div style={{ fontSize: '10px', color: '#6b7280', textTransform: 'uppercase', letterSpacing: '0.05em' }}>
          {data.statusText}
        </div>
        <div style={{ textAlign: 'right' }}>
          <div style={{ fontSize: '16px', fontWeight: 'bold', color: '#111827', lineHeight: '1' }}>
            {data.count}
          </div>
          {data.removed > 0 && (
            <div style={{ fontSize: '11px', color: '#ef4444', marginTop: '2px' }}>
              -{data.removed}
            </div>
          )}
        </div>
      </div>
      
      <Handle type="source" position={Position.Bottom} style={{ background: '#555' }} />
    </div>
  );
};

const nodeTypes = { custom: PipelineNode };

const CurationDag = (props: ComponentProps) => {
  const { args } = props;
  const payload = args.payload || { nodes: [], edges: [] };
  
  const nodes = useMemo(() => {
    return payload.nodes.map((n: any) => ({
      id: n.id,
      type: 'custom',
      position: { x: n.x * 240, y: n.y * 120 },
      data: {
        label: n.label,
        glyph: n.glyph,
        statusText: n.statusText,
        color: n.color,
        count: n.count,
        removed: n.removed,
      },
    }));
  }, [payload.nodes]);

  const edges = useMemo(() => {
    return payload.edges.map((e: any, idx: number) => ({
      id: `e-${e.source}-${e.target}-${idx}`,
      source: e.source,
      target: e.target,
      type: 'smoothstep',
      animated: true,
      style: { stroke: '#9ca3af', strokeWidth: 2 },
      markerEnd: {
        type: MarkerType.ArrowClosed,
        color: '#9ca3af',
      },
    }));
  }, [payload.edges]);

  const onNodeClick = (_: React.MouseEvent, node: any) => {
    Streamlit.setComponentValue(node.id);
  };

  useEffect(() => {
    Streamlit.setFrameHeight(600);
  }, []);

  return (
    <div style={{ width: '100%', height: '600px', border: '1px solid #e5e7eb', borderRadius: '12px', background: '#f9fafb' }}>
      <ReactFlow 
        nodes={nodes} 
        edges={edges} 
        nodeTypes={nodeTypes}
        onNodeClick={onNodeClick}
        fitView
        fitViewOptions={{ padding: 0.2 }}
        nodesConnectable={false}
        nodesDraggable={true}
        elementsSelectable={true}
      >
        <Background color="#d1d5db" gap={20} />
        <Controls showInteractive={false} />
      </ReactFlow>
    </div>
  );
};

export default withStreamlitConnection(CurationDag);
